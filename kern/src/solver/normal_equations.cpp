#include "normal_equations.h"

#include "run/progress.h"
#include "util/introsort.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <utility>

namespace areko::remesher::solver {

namespace {

constexpr double keepThreshold = 1e-16;

} // namespace

NormalEquations::Result NormalEquations::assemble(const ListMatrix& rows, std::span<const double> weights,
                                                  std::span<const double> targets, Scope scope,
                                                  const run::Progress& progress)
{
    if (scope != Scope::RightHandSideOnly) {
        if (scope == Scope::All) {
            transposed_ = CompactMatrix::fromListsTransposed(rows);
            progress.check();
        }
        buildLower(rows, weights, progress);
        if (matrix_.nonZeros() == 0) {
            return Result::EmptyMatrix;
        }
        progress.check();
        mirrorToFull();
        progress.check();
    }
    buildRightHandSide(weights, targets);
    return Result::Ok;
}

void NormalEquations::buildLower(const ListMatrix& rows, std::span<const double> weights,
                                 const run::Progress& progress)
{
    const int n = rows.otherDimension();
    std::vector<double> sum(static_cast<std::size_t>(n), 0.0);
    std::vector<std::uint8_t> touched(static_cast<std::size_t>(n), 0);
    std::vector<int> order;
    matrix_ = CompactMatrix(n);
    for (int c = 0; c < n; ++c) {
        progress.check();
        order.clear();
        for (int k = transposed_.columnStart(c); k < transposed_.columnEnd(c); ++k) {
            const int r = transposed_.row(k);
            const double arc = transposed_.value(k);
            const double wr = weights[static_cast<std::size_t>(r)];
            // Row r from its end while the unknown is not left of c (rows are sorted).
            const ListMatrix::List& row = rows.list(r);
            for (auto it = row.rbegin(); it != row.rend() && it->index >= c; ++it) {
                if (it->value == 0.0) {
                    continue;
                }
                const auto j = static_cast<std::size_t>(it->index);
                if (touched[j] == 0) {
                    touched[j] = 1;
                    order.push_back(it->index);
                }
                sum[j] = it->value * (wr * arc) + sum[j];
            }
        }
        util::Introsort::sort(order.begin(), order.end(), [](int a, int b) { return a < b; });
        for (const int j : order) {
            const auto u = static_cast<std::size_t>(j);
            if (std::fabs(sum[u]) > keepThreshold) {
                matrix_.append(j, sum[u]);
            }
            sum[u] = 0.0;
            touched[u] = 0;
        }
        matrix_.closeColumn();
    }
}

void NormalEquations::mirrorToFull()
{
    // Column c of the transposed lower triangle holds the rows j <= c ascending,
    // the diagonal last when it is stored.
    const CompactMatrix upper = matrix_.transposed();
    CompactMatrix full(matrix_.rows());
    for (int c = 0; c < matrix_.columns(); ++c) {
        const int begin = matrix_.columnStart(c);
        const int end = matrix_.columnEnd(c);
        const bool diagonalFirst = begin < end && matrix_.row(begin) == c;
        const int upperEnd = upper.columnEnd(c) - (diagonalFirst ? 1 : 0);
        for (int k = upper.columnStart(c); k < upperEnd; ++k) {
            full.append(upper.row(k), upper.value(k));
        }
        for (int k = begin; k < end; ++k) {
            full.append(matrix_.row(k), matrix_.value(k));
        }
        full.closeColumn();
    }
    matrix_ = std::move(full);
}

void NormalEquations::buildRightHandSide(std::span<const double> weights, std::span<const double> targets)
{
    const int n = transposed_.columns();
    rhs_.assign(static_cast<std::size_t>(n), 0.0);
    for (int c = 0; c < n; ++c) {
        double g = 0.0;
        for (int k = transposed_.columnStart(c); k < transposed_.columnEnd(c); ++k) {
            const auto r = static_cast<std::size_t>(transposed_.row(k));
            g = g + (weights[r] * targets[r]) * transposed_.value(k);
        }
        rhs_[static_cast<std::size_t>(c)] = g;
    }
}

NormalEquations::ListSystem NormalEquations::assembleLists(const ListMatrix& rows, std::span<const double> weights,
                                                           std::span<const double> targets,
                                                           const run::Progress& progress)
{
    const int n = rows.otherDimension();
    // Column c of A: rows ascending, each in storage order; column r of its
    // transpose: the unknowns of row r, ascending.
    const CompactMatrix columns = CompactMatrix::fromListsTransposed(rows);
    const CompactMatrix rowUnknowns = columns.transposed();

    ListSystem out{ListMatrix(n, n), std::vector<double>(static_cast<std::size_t>(n), 0.0)};
    std::vector<int> pattern;
    std::vector<int> merged;
    for (int c = 0; c < n; ++c) {
        pattern.clear();
        for (int k = columns.columnStart(c); k < columns.columnEnd(c); ++k) {
            const std::span<const int> add = rowUnknowns.columnRows(columns.row(k));
            merged.clear();
            std::set_union(pattern.begin(), pattern.end(), add.begin(), add.end(), std::back_inserter(merged));
            pattern.swap(merged);
        }
        ListMatrix::List& target = out.matrix.list(c);
        target.reserve(pattern.size());
        for (const int j : pattern) {
            double v = 0.0;
            int p = columns.columnStart(c);
            int q = columns.columnStart(j);
            while (p < columns.columnEnd(c) && q < columns.columnEnd(j)) {
                if (columns.row(p) < columns.row(q)) {
                    ++p;
                } else if (columns.row(q) < columns.row(p)) {
                    ++q;
                } else {
                    v = v + (columns.value(p) * columns.value(q)) * weights[static_cast<std::size_t>(columns.row(p))];
                    ++p;
                    ++q;
                }
            }
            target.push_back({j, v});
        }
        progress.check();
    }
    for (int j = 0; j < n; ++j) {
        double g = 0.0;
        for (int k = columns.columnStart(j); k < columns.columnEnd(j); ++k) {
            const auto i = static_cast<std::size_t>(columns.row(k));
            g = g + (columns.value(k) * targets[i]) * weights[i];
        }
        out.rhs[static_cast<std::size_t>(j)] = g;
    }
    return out;
}

} // namespace areko::remesher::solver
