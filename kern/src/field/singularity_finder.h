#pragma once

#include "cross_field.h"

#include "mesh/mesh.h"
#include "prep/symmetry_cache.h"

#include <vector>

namespace areko::remesher::field {

// Finds the points around which the crosses do not close up.
class SingularityFinder {
public:
    // Walks the triangle fan of every inner point: the quarter turns across
    // its edges add up to the index, the angles of the field across them to
    // the turn of the field, the corner angles to the turn of the surface.
    // A point is listed when its index is not 0 modulo 4, or when field and
    // surface turn alike (one whole turn of the crosses). Border points and
    // points without triangles are left out.
    // `mirror` are the pairs of a symmetric part (or null); with them and a
    // valid evaluation (`mirrorValid`), a whole turn is not listed at a point
    // of side 0 that has a mirror point other than itself.
    [[nodiscard]] static Singularities find(const mesh::Mesh& mesh, const CrossField& field,
                                            const std::vector<float>& connection,
                                            const prep::SymmetryCache* mirror, bool mirrorValid);
};

} // namespace areko::remesher::field
