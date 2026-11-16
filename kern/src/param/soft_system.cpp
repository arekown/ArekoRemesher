#include "soft_system.h"

#include "param/elimination.h"
#include "run/progress.h"

#include <cstddef>
#include <utility>

namespace areko::remesher::param {

namespace {

// Weight of the rows that tie every unknown to its present value.
constexpr double tieWeight = static_cast<double>(1e-4f);
static_assert(tieWeight == 9.9999997473787516e-05);

} // namespace

void SoftSystem::addSoft(TermList row, double target, double weight)
{
    soft_.push_back(std::move(row));
    softTargets_.push_back(target);
    softWeights_.push_back(weight);
}

void SoftSystem::addHard(TermList row, double target)
{
    hard_.push_back(std::move(row));
    hardTargets_.push_back(target);
}

ReducedSoftSystem SoftSystem::reduce(const System& system, std::span<const int> preferred, std::span<const double> x,
                                     const run::Progress& progress) &&
{
    // The conditions of the main system leave every row in its normal form.
    system.constraints.substituteInto(soft_);
    progress.check();
    system.constraints.substituteInto(hard_);
    // A hard row that lost all its terms goes, with its target.
    std::size_t kept = 0;
    for (std::size_t i = 0; i < hard_.size(); ++i) {
        if (hard_[i].empty()) {
            continue;
        }
        if (kept != i) {
            hard_[kept] = std::move(hard_[i]);
            hardTargets_[kept] = hardTargets_[i];
        }
        ++kept;
    }
    hard_.resize(kept);
    hardTargets_.resize(kept);
    // The targets become terms of the constant column.
    bool constants = false;
    for (std::size_t i = 0; i < hard_.size(); ++i) {
        if (hardTargets_[i] != 0.0) {
            hard_[i].push_back({constantColumn(), hardTargets_[i]});
            hardTargets_[i] = 0.0;
            constants = true;
        }
    }

    EliminationSetup setup;
    setup.mode = EliminationMode::SoftSystem;
    setup.preferred = preferred;
    setup.constantColumn = constantColumn();
    Elimination elimination(hard_, columnCount_, setup, progress);
    (void)elimination.run();
    progress.check();

    ReducedSoftSystem out;
    out.conditions = ConstraintSet(columnCount_);
    elimination.enterPivots(out.conditions.variableMap(), 0);
    out.conditions.append(elimination.takeConstraints());
    out.conditions.substituteInto(soft_);
    if (constants) {
        // What a soft row holds of the constant column moves to its target.
        for (std::size_t i = 0; i < soft_.size(); ++i) {
            TermList& row = soft_[i];
            for (std::size_t k = 0; k < row.size();) {
                if (row[k].index == constantColumn()) {
                    softTargets_[i] = softTargets_[i] + row[k].value;
                    row.erase(row.begin() + static_cast<std::ptrdiff_t>(k));
                } else {
                    ++k;
                }
            }
        }
    }

    // The columns left in the soft rows become the unknowns.
    out.unknownOf.assign(static_cast<std::size_t>(columnCount_), -1);
    for (TermList& row : soft_) {
        for (Term& term : row) {
            int& unknown = out.unknownOf[static_cast<std::size_t>(term.index)];
            if (unknown == -1) {
                unknown = static_cast<int>(out.columnOf.size());
                out.columnOf.push_back(term.index);
            }
            term.index = unknown;
        }
    }
    for (int column = 0; column < columnCount_; ++column) {
        const int unknown = out.unknownOf[static_cast<std::size_t>(column)];
        if (unknown == -1) {
            continue;
        }
        double target = 1.0;
        if (column < variableCount_) {
            const int mainUnknown = system.variableToUnknown[static_cast<std::size_t>(column)];
            target = mainUnknown < 0 ? 0.0 : x[static_cast<std::size_t>(mainUnknown)];
        }
        addSoft(TermList{{unknown, 1.0}}, target, tieWeight);
    }

    out.rows = solver::ListMatrix(static_cast<int>(soft_.size()), out.unknownCount());
    for (std::size_t i = 0; i < soft_.size(); ++i) {
        out.rows.list(static_cast<int>(i)) = std::move(soft_[i]);
    }
    out.weights = std::move(softWeights_);
    out.targets = std::move(softTargets_);
    return out;
}

} // namespace areko::remesher::param
