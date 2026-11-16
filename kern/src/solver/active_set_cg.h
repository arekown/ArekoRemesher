#pragma once

#include "solver/sparse_matrix.h"

#include <span>
#include <vector>

namespace areko::remesher::solver {

// Conjugate gradients without preconditioning on a symmetric list system,
// starting from x and touching only an active set that grows by one ring of
// neighbours per iteration (unknowns outside it keep their value).
class ActiveSetCg {
public:
    struct Options {
        double tolerance = 0.1;   // relative to the initial residual
        int maxIterations = 100;  // -1 = 1000
    };

    // `set` is the initial active set (order matters).
    // Returns true when the residual over the active set fell below the
    // tolerance (also when it already was, or the set is empty).
    [[nodiscard]] static bool solve(const Options& options, std::vector<double>& x, std::vector<int> set,
                                    const ListMatrix& matrix, std::span<const double> rhs);
};

} // namespace areko::remesher::solver
