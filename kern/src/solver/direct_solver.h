#pragma once

#include "solver/sparse_matrix.h"

#include <span>
#include <vector>

namespace areko::remesher::run {
class Progress;
}

namespace areko::remesher::solver {

// Direct solver for symmetric sparse systems M x = b (full storage, only the
// lower triangle is read): LDLᵀ with an approximate minimum degree ordering.
// Every call orders, analyses and factors anew; nothing is reused.
class DirectSolver {
public:
    // Compressed input. `dropZeros`: entries with value 0.0 are left out of the
    // pattern (NaN stays), otherwise every stored entry counts.
    bool solve(const CompactMatrix& matrix, bool dropZeros, std::span<const double> rhs,
               const run::Progress& progress, std::vector<double>& x);
    // List input: list c is column c, every entry counts (also 0.0).
    bool solve(const ListMatrix& matrix, std::span<const double> rhs, const run::Progress& progress,
               std::vector<double>& x);

    // Elimination order of the last factorization: position k holds the
    // unknown eliminated k-th (the inverse permutation of the ordering).
    [[nodiscard]] const std::vector<int>& ordering() const noexcept { return ordering_; }

private:
    // `matrix` is square; every stored entry is part of the pattern.
    bool factorAndSolve(const CompactMatrix& matrix, std::span<const double> rhs, const run::Progress& progress,
                        std::vector<double>& x);

    std::vector<int> ordering_;
};

} // namespace areko::remesher::solver
