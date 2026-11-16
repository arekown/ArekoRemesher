#include "constrained_driver.h"

#include "run/progress.h"
#include "util/rounding.h"

#include <cfloat>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace areko::remesher::solver {

namespace {

// A solution is usable when it is not empty and every value is finite.
bool usable(const std::vector<double>& x)
{
    if (x.empty()) {
        return false;
    }
    for (const double v : x) {
        if (!std::isfinite(v)) {
            return false;
        }
    }
    return true;
}

// Stage tolerances of the batch rounding, as fractions of the budget.
constexpr double batchStages[] = {1e-5, 1e-4, 1e-3, 1e-2, 0.05};

} // namespace

ConstrainedDriver::ConstrainedDriver(CompactMatrix& matrix, std::vector<double>& rhs, DriverCallback& callback,
                                     const run::Progress& progress, const DriverOptions& options)
    : matrix_(&matrix), rhs_(rhs), callback_(callback), progress_(progress), options_(options)
{
    if (options_.localUpdate) {
        throw std::invalid_argument("local updates need a list system");
    }
}

ConstrainedDriver::ConstrainedDriver(ListMatrix& matrix, std::vector<double>& rhs, DriverCallback& callback,
                                     const run::Progress& progress, const DriverOptions& options)
    : matrix_(&matrix), rhs_(rhs), callback_(callback), progress_(progress), options_(options)
{
}

void ConstrainedDriver::setIntegerUnknowns(std::vector<int> unknowns, std::vector<double> periods)
{
    unknowns_ = std::move(unknowns);
    periods_ = std::move(periods);
}

bool ConstrainedDriver::solveSystem(std::vector<double>& x)
{
    return std::visit([this, &x](const auto* matrix) { return solveIn(*matrix, x); }, matrix_);
}

bool ConstrainedDriver::solveIn(const CompactMatrix& matrix, std::vector<double>& x)
{
    return direct_.solve(matrix, zeroCounter_ >= 1, rhs_, progress_, x);
}

bool ConstrainedDriver::solveIn(const ListMatrix& matrix, std::vector<double>& x)
{
    return direct_.solve(matrix, rhs_, progress_, x);
}

bool ConstrainedDriver::run(std::vector<double>& x)
{
    x.clear();
    status_ = Status::Ok;
    zeroCounter_ = -1;
    const bool solved = solveSystem(x);
    zeroCounter_ = 0;
    if (!solved || !usable(x)) {
        status_ = Status::SolveFailed;
        return false;
    }
    progress_.check();
    if (options_.localUpdate) {
        queue_.reset(rhs_.size());
    }
    const int n = static_cast<int>(unknowns_.size());
    fixed_ = mesh::BitArray(static_cast<std::size_t>(n));
    int fixedCount = 0;
    int changedRounds = 0;
    int round = -1;
    int sinceFullSolve = 0;
    batchMemo_ = 0.0;
    DriverRequest request;
    request.x = &x;
    std::vector<double> next;
    while (true) {
        if (fixedCount >= n && !options_.keepGoing) {
            return true;
        }
        if (request.finished) {
            return true;
        }
        progress_.check();
        ++round;
        request.extra = nullptr;
        request.answer = DriverRequest::Answer::None;
        request.systemChanged = false;
        request.selection.clear();
        request.allowLocal = true;
        request.noMark = false;
        request.changedRounds = changedRounds;
        request.round = round;
        request.fixedCount = fixedCount;
        callback_.onRequest(request);

        if (request.answer == DriverRequest::Answer::Stop) {
            return false;
        }
        const bool delivered = request.answer == DriverRequest::Answer::Delivered;
        if (delivered && !options_.keepGoing && request.selection.empty()) {
            if (!request.systemChanged) {
                return false;
            }
            ++changedRounds;
        } else if (options_.keepGoing && !request.finished) {
            if (!request.systemChanged || request.extra == nullptr || request.extra->empty()) {
                return false;
            }
            ++changedRounds;
        } else if (request.systemChanged) {
            ++changedRounds;
        }

        if (delivered) {
            if (!request.noMark) {
                for (const DriverRequest::Selection& chosen : request.selection) {
                    const auto p = static_cast<std::size_t>(chosen.position);
                    if (p >= fixed_.size()) {
                        fixed_.resize(p + 1);
                    }
                    fixed_.set(p);
                }
            }
        } else {
            progress_.check();
            request.selection.clear();
            if (options_.batchMode) {
                pickBatch(x, request);
                if (request.selection.empty()) {
                    return false;
                }
            } else if (!pickGreedy(x, request)) {
                return false;
            }
        }
        if (!request.selection.empty()) {
            progress_.check();
            callback_.onSelection(request);
        }
        progress_.check();
        for (std::size_t k = 0; k < request.selection.size(); ++k) {
            const int var = unknowns_[static_cast<std::size_t>(request.selection[k].position)];
            const double value = request.selection[k].value;
            ++fixedCount;
            x[static_cast<std::size_t>(var)] = value;
            if (options_.localUpdate) {
                for (const SparseEntry& e : std::get<ListMatrix*>(matrix_)->list(var)) {
                    if (e.value != 0.0 && e.index != var) {
                        queue_.pushBack(e.index);
                    }
                }
            }
            fixVariable(var, value);
            if (k + 1 < request.selection.size()) {
                progress_.check();
            }
        }
        progress_.check();

        bool local = false;
        if (options_.localUpdate) {
            local = true;
            if (options_.fullSolveEvery > 0) {
                ++sinceFullSolve;
                local = sinceFullSolve % options_.fullSolveEvery != 0;
            }
            if (!request.allowLocal) {
                local = false;
            }
        }
        if (!options_.keepGoing && fixedCount == n) {
            callback_.onAllFixed();
            queue_.clear();
        } else if (!local) {
            queue_.clear();
        } else if (localUpdate(request, x)) {
            continue;
        }
        progress_.check();
        if (!solveSystem(next) || !usable(next)) {
            status_ = Status::SolveFailed;
            return false;
        }
        x.swap(next);
        sinceFullSolve = 0;
    }
}

bool ConstrainedDriver::pickGreedy(const std::vector<double>& x, DriverRequest& request)
{
    // Greedy rounding goes to the nearest integer; periods do not count.
    double error = 0.0;
    double value = 0.0;
    const int pick = bestFree(x, false, error, value);
    if (pick == -1) {
        return false;
    }
    take(request, pick, value);
    return true;
}

double ConstrainedDriver::nearestWhole(const std::vector<double>& x, int position) const
{
    const double v = x[static_cast<std::size_t>(unknowns_[static_cast<std::size_t>(position)])];
    return static_cast<double>(util::roundHalfAway(v));
}

void ConstrainedDriver::roundingOf(const std::vector<double>& x, int position, bool withPeriods, double& error,
                                   double& value) const
{
    const auto p = static_cast<std::size_t>(position);
    const double v = x[static_cast<std::size_t>(unknowns_[p])];
    if (!withPeriods || periods_.empty()) {
        value = nearestWhole(x, position);
        error = std::fabs(v - value);
        return;
    }
    const PeriodRounding rounded = toPeriod(v, periods_[p]);
    error = rounded.error;
    value = rounded.value;
}

ConstrainedDriver::PeriodRounding ConstrainedDriver::toPeriod(double value, double period) noexcept
{
    const double y = value / period;
    const int r = util::roundHalfAway(y);
    return {period * static_cast<double>(r), std::fabs(y - static_cast<double>(r))};
}

int ConstrainedDriver::bestFree(const std::vector<double>& x, bool withPeriods, double& error, double& value) const
{
    int pick = -1;
    error = DBL_MAX;
    value = 0.0;
    const int n = static_cast<int>(unknowns_.size());
    for (int p = 0; p < n; ++p) {
        if (!isFree(p)) {
            continue;
        }
        double e = 0.0;
        double v = 0.0;
        roundingOf(x, p, withPeriods, e, v);
        if (e < error) {
            error = e;
            value = v;
            pick = p;
        }
    }
    return pick;
}

void ConstrainedDriver::take(DriverRequest& request, int position, double value)
{
    request.selection.push_back({position, value});
    fixed_.set(static_cast<std::size_t>(position));
}

void ConstrainedDriver::batchStage(const std::vector<double>& x, double tolerance, double& batchError,
                                   DriverRequest& request)
{
    const int n = static_cast<int>(unknowns_.size());
    for (int p = 0; p < n; ++p) {
        if (!isFree(p)) {
            continue;
        }
        double error = 0.0;
        double value = 0.0;
        roundingOf(x, p, true, error, value);
        if (error < tolerance) {
            batchError = error + batchError;
            if (batchError > options_.batchBudget) {
                return;
            }
            take(request, p, value);
        }
    }
}

void ConstrainedDriver::pickBatch(const std::vector<double>& x, DriverRequest& request)
{
    const double budget = options_.batchBudget;
    double batchError = 0.0;
    if (batchMemo_ == 0.0) {
        // Stages with growing tolerance while the batch stays below the budget.
        for (const double fraction : batchStages) {
            if (fraction != batchStages[0] && !(batchError < budget)) {
                break;
            }
            batchStage(x, budget * fraction, batchError, request);
            batchMemo_ = budget * fraction;
        }
    }
    const int n = static_cast<int>(unknowns_.size());
    while (static_cast<int>(request.selection.size()) < n) {
        double best = 0.0;
        double bestValue = 0.0;
        const int pick = bestFree(x, true, best, bestValue);
        if (pick == -1) {
            return;
        }
        batchError = batchError + best;
        if (budget < batchError) {
            // One position is always taken.
            if (request.selection.empty()) {
                take(request, pick, bestValue);
            }
            batchMemo_ = best;
            return;
        }
        take(request, pick, bestValue);
    }
}

void ConstrainedDriver::fixVariable(int var, double value)
{
    std::visit([this, var, value](auto* matrix) { fixIn(*matrix, var, value); }, matrix_);
}

void ConstrainedDriver::fixIn(CompactMatrix& m, int var, double value)
{
    const auto j = static_cast<std::size_t>(var);
    for (int k = m.columnStart(var); k < m.columnEnd(var); ++k) {
        const int i = m.row(k);
        if (i == var) {
            m.value(k) = 1.0;
            rhs_[j] = value;
            continue;
        }
        if (m.value(k) != 0.0) {
            ++zeroCounter_;
            m.value(k) = 0.0;
        }
        // Mirror entry (row var in column i); always present in a symmetric pattern.
        int mirror = -1;
        for (int q = m.columnStart(i); q < m.columnEnd(i); ++q) {
            if (m.row(q) == var) {
                mirror = q;
                break;
            }
        }
        if (mirror < 0) {
            continue;
        }
        const double a = m.value(mirror);
        if (a != 0.0) {
            ++zeroCounter_;
            m.value(mirror) = 0.0;
        }
        rhs_[static_cast<std::size_t>(i)] = rhs_[static_cast<std::size_t>(i)] - a * value;
    }
}

void ConstrainedDriver::fixIn(ListMatrix& m, int var, double value)
{
    for (const SparseEntry& e : m.list(var)) {
        if (e.index == var) {
            continue;
        }
        ListMatrix::List& other = m.list(e.index);
        const int at = ListMatrix::find(other, var);
        if (at < 0) {
            continue;
        }
        const double a = other[static_cast<std::size_t>(at)].value;
        other.erase(other.begin() + at);
        rhs_[static_cast<std::size_t>(e.index)] = rhs_[static_cast<std::size_t>(e.index)] - a * value;
    }
    m.list(var).assign(1, SparseEntry{var, 1.0});
    rhs_[static_cast<std::size_t>(var)] = value;
}

bool ConstrainedDriver::localUpdate(const DriverRequest& request, std::vector<double>& x)
{
    if (options_.keepGoing && request.extra != nullptr) {
        for (const int i : *request.extra) {
            queue_.pushFront(i);
        }
    }
    std::vector<int> set;
    queue_.drainTo(set);
    return ActiveSetCg::solve(options_.local, x, std::move(set), *std::get<ListMatrix*>(matrix_), rhs_);
}

} // namespace areko::remesher::solver
