#pragma once

#include "solver/sparse_matrix.h"

#include <span>
#include <vector>

namespace areko::remesher::run {
class Progress;
}

namespace areko::remesher::solver {

// Normal equations H x = g with H = Aᵀ W A and g = Aᵀ W b of a weighted least
// squares system given as row lists A (index = unknown, every row sorted by
// index, otherDimension = number of unknowns), weights w and targets b.
class NormalEquations {
public:
    enum class Scope {
        All,                // Aᵀ, H and g
        KeepTranspose,      // H and g from the kept Aᵀ (new weights, same rows)
        RightHandSideOnly,  // g only
    };
    enum class Result { Ok, EmptyMatrix };

    // Compressed path: H keeps only entries with |h| > 1e-16, is built as lower
    // triangle and then mirrored to full storage, every column ascending.
    // EmptyMatrix: H has no entry; g is not computed then.
    Result assemble(const ListMatrix& rows, std::span<const double> weights, std::span<const double> targets,
                    Scope scope, const run::Progress& progress);

    [[nodiscard]] const CompactMatrix& transposed() const noexcept { return transposed_; }
    [[nodiscard]] CompactMatrix& matrix() noexcept { return matrix_; }
    [[nodiscard]] const CompactMatrix& matrix() const noexcept { return matrix_; }
    [[nodiscard]] std::vector<double>& rightHandSide() noexcept { return rhs_; }
    [[nodiscard]] const std::vector<double>& rightHandSide() const noexcept { return rhs_; }

    // List path: H as one sorted list per unknown in full storage, without a
    // threshold (entries that sum to zero stay), with its right-hand side.
    struct ListSystem {
        ListMatrix matrix;
        std::vector<double> rhs;
    };
    [[nodiscard]] static ListSystem assembleLists(const ListMatrix& rows, std::span<const double> weights,
                                                  std::span<const double> targets, const run::Progress& progress);

private:
    void buildLower(const ListMatrix& rows, std::span<const double> weights, const run::Progress& progress);
    void mirrorToFull();
    void buildRightHandSide(std::span<const double> weights, std::span<const double> targets);

    CompactMatrix transposed_;
    CompactMatrix matrix_;
    std::vector<double> rhs_;
};

} // namespace areko::remesher::solver
