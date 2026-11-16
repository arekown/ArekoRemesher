#pragma once

#include "curves/guide_curve.h"
#include "mesh/bit_array.h"
#include "mesh/mesh.h"

#include <span>
#include <vector>

namespace areko::remesher::field {

// Curvature carried by the edges of a triangle mesh: per edge its signed
// dihedral angle, its unit direction and the tensor length * angle * e e^T.
class EdgeCurvature {
public:
    // Per edge with end points a, b and length l:
    //  - direction (P[b] - P[a]) / l (zero for l = 0),
    //  - angle 0 unless the edge has exactly two corners; then, with the face
    //    normals n0, n1 of its first and second corner (double) and their
    //    opposite points o0, o1: s = (P[o1] - P[o0]) . n0 and
    //    cos = n0 . n1 / (|n0| |n1|) in double; s < 0 (convex) gives
    //    pi / 0 / acos(cos) for cos <= -1 / >= 1 / else, otherwise the
    //    negative of that (-0 for >= 1); a zero denominator gives 0 resp. -0,
    //    a NaN cosine stays NaN,
    //  - tensor [xx, xy, xz, yy, yz, zz] of the direction times
    //    w = double(l * angle).
    explicit EdgeCurvature(const mesh::Mesh& mesh);

    [[nodiscard]] const std::vector<float>& lengths() const noexcept { return lengths_; }
    [[nodiscard]] const std::vector<float>& angles() const noexcept { return angles_; }
    [[nodiscard]] const std::vector<mesh::Vec3f>& directions() const noexcept { return directions_; }
    [[nodiscard]] std::span<const double, 6> tensor(int e) const
    {
        return std::span<const double, 6>(tensors_.data() + 6 * static_cast<std::size_t>(e), 6);
    }

    // End points of the edges between successive points of the cut chains
    // of active crease curves whose angle is not below 20 degrees (NaN counts).
    [[nodiscard]] mesh::BitArray sharpPoints(const mesh::Mesh& mesh,
                                             const std::vector<curves::GuideCurve>& curves) const;

private:
    std::vector<float> lengths_;
    std::vector<float> angles_;
    std::vector<mesh::Vec3f> directions_;
    std::vector<double> tensors_;
};

} // namespace areko::remesher::field
