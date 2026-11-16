#pragma once

#include "mesh/mesh.h"

#include <span>

namespace areko::remesher::prep {

// Point and face pairs of a mirror-symmetric mesh (-1 = no partner).
struct MirrorLists {
    std::span<const int> points;
    std::span<const int> faces;
};

// Turns a polygon or quad mesh into a triangle mesh.
//  - quads: a corner whose point has exactly two corners and two edges starts
//    the split; otherwise, with `angleRule`, the corner angle sums and a parity
//    vote over the points decide the diagonal, without it the shorter one,
//  - other polygons: PolygonTriangulator,
//  - with `mirror`, faces follow the start point of their mirror face (no
//    diagonal rule, sharp edges are dropped); with `mirror` but no face
//    pairs, quads take the shorter diagonal (no angle rule).
// Points and point maps move over, face maps repeat per triangle, sharp edges
// keep their polygon edges (new diagonals are never sharp), UV sets follow the
// corners. Other maps and given normals are not taken over. A triangle mesh
// comes back unchanged.
class Triangulator {
public:
    struct Options {
        bool angleRule = true;
        const MirrorLists* mirror = nullptr;
    };

    [[nodiscard]] static mesh::Mesh run(mesh::Mesh&& source, const Options& options);
};

} // namespace areko::remesher::prep
