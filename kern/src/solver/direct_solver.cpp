#include "direct_solver.h"

#include "run/progress.h"

#include <Eigen/OrderingMethods>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include <algorithm>

namespace areko::remesher::solver {

namespace {

using EigenMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
// Lower triangle as input; AMD on the pattern of the full symmetric matrix;
// solving multiplies by the reciprocal of D.
using EigenLdlt = Eigen::SimplicialLDLT<EigenMatrix, Eigen::Lower, Eigen::AMDOrdering<int>>;

} // namespace

bool DirectSolver::solve(const CompactMatrix& matrix, bool dropZeros, std::span<const double> rhs,
                         const run::Progress& progress, std::vector<double>& x)
{
    if (matrix.rows() != matrix.columns()) {
        return false;
    }
    if (!dropZeros) {
        return factorAndSolve(matrix, rhs, progress, x);
    }
    CompactMatrix kept(matrix.rows());
    for (int c = 0; c < matrix.columns(); ++c) {
        for (int k = matrix.columnStart(c); k < matrix.columnEnd(c); ++k) {
            if (matrix.value(k) != 0.0) {
                kept.append(matrix.row(k), matrix.value(k));
            }
        }
        kept.closeColumn();
    }
    return factorAndSolve(kept, rhs, progress, x);
}

bool DirectSolver::solve(const ListMatrix& matrix, std::span<const double> rhs, const run::Progress& progress,
                         std::vector<double>& x)
{
    if (!matrix.square()) {
        return false;
    }
    CompactMatrix columns(matrix.listCount());
    for (int c = 0; c < matrix.listCount(); ++c) {
        for (const SparseEntry& e : matrix.list(c)) {
            columns.append(e.index, e.value);
        }
        columns.closeColumn();
    }
    return factorAndSolve(columns, rhs, progress, x);
}

bool DirectSolver::factorAndSolve(const CompactMatrix& matrix, std::span<const double> rhs,
                                  const run::Progress& progress, std::vector<double>& x)
{
    const int n = matrix.columns();
    if (rhs.size() != static_cast<std::size_t>(n)) {
        return false;
    }
    if (n == 0) {
        ordering_.clear();
        x.clear();
        return true;
    }
    EigenMatrix input(n, n);
    input.resizeNonZeros(static_cast<Eigen::Index>(matrix.nonZeros()));
    std::copy(matrix.starts().begin(), matrix.starts().end(), input.outerIndexPtr());
    std::copy(matrix.rowIndices().begin(), matrix.rowIndices().end(), input.innerIndexPtr());
    std::copy(matrix.values().begin(), matrix.values().end(), input.valuePtr());
    progress.check();

    EigenLdlt ldlt;
    ldlt.analyzePattern(input);
    if (ldlt.info() != Eigen::Success) {
        return false;
    }
    progress.check();
    ldlt.factorize(input);
    if (ldlt.info() != Eigen::Success) {
        return false;
    }
    progress.check();
    const Eigen::Map<const Eigen::VectorXd> b(rhs.data(), n);
    const Eigen::VectorXd z = ldlt.solve(b);
    if (ldlt.info() != Eigen::Success) {
        return false;
    }
    const auto& order = ldlt.permutationPinv().indices();
    ordering_.assign(order.data(), order.data() + order.size());
    x.assign(z.data(), z.data() + z.size());
    return true;
}

} // namespace areko::remesher::solver
