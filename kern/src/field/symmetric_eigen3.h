#pragma once

#include "mesh/vec3.h"

#include <array>

namespace areko::remesher::field {

// Eigen decomposition of a symmetric 3x3 matrix, all in double.
class SymmetricEigen3 {
public:
    struct Pair {
        mesh::Vec3d vector;
        double value = 0.0;
    };

    // Matrix as [xx, xy, xz, yy, yz, zz]. Cyclic Jacobi with threshold: stops
    // when the off-diagonal sum is below 1e-20 or after 70 sweeps; the
    // threshold is sum * 0.2 / 9 in the first three sweeps, 0 after; from the
    // fifth sweep on an entry whose hundredfold is below 1e-20 is set to 0.
    // The pairs come sorted by |value| ascending (compare and swap (0,1),
    // (1,2), (0,1), only on a strictly smaller value), so equal values keep
    // the order xx, yy, zz; signs of the vectors are left as they come.
    [[nodiscard]] static std::array<Pair, 3> solve(const std::array<double, 6>& m);
};

} // namespace areko::remesher::field
