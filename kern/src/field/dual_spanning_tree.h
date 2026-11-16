#pragma once

#include "mesh/bit_array.h"
#include "mesh/mesh.h"

namespace areko::remesher::field {

// Spanning forest over the triangles (neighbours across inner edges), grown
// from the anchors along the shortest paths between triangle centres.
class DualSpanningTree {
public:
    // One bit per edge: the free edges, i.e. the inner edges that are neither
    // a tree edge nor lie between two anchors. Along a free edge the quarter
    // turn between the two crosses is still to be chosen.
    // Anchors in `mirrorHalf` do not start a tree; without a start triangle 0 does.
    [[nodiscard]] static mesh::BitArray freeEdges(const mesh::Mesh& mesh, const mesh::BitArray& anchors,
                                                  const mesh::BitArray& mirrorHalf);
};

} // namespace areko::remesher::field
