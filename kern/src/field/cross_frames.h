#pragma once

#include "cross_field.h"

#include "mesh/mesh.h"
#include "prep/symmetry_cache.h"

#include <vector>

namespace areko::remesher::field {

// Frames of the triangles and the angle that carries one frame into the next.
class CrossFrames {
public:
    // Frame of every triangle: the first axis runs along its reference edge
    // (corner 0 to corner 1), the second is normal x first axis. With the
    // mirror pairs of a symmetric part, a triangle on side 0 that has a mirror
    // triangle takes the mirror image of its partner's reference edge instead,
    // so both frames mirror onto each other. Fills frame, start corner and
    // reference length of `field`.
    static void compute(const mesh::Mesh& mesh, prep::SymmetryCache* mirror, CrossField& field);

    // Per inner edge the angle (-pi, pi] under which the first axis of its
    // second triangle appears in the frame of its first triangle, after
    // turning the second triangle into the plane of the first about their
    // common edge; counted about the normal of the first. Other edges get 0.
    [[nodiscard]] static std::vector<float> connection(const mesh::Mesh& mesh,
                                                       const std::vector<mesh::Vec3f>& frameE1);
};

} // namespace areko::remesher::field
