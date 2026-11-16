#pragma once

#include "field/cut_copy.h"
#include "mesh/adjacency.h"

#include <optional>
#include <span>
#include <vector>

namespace areko::remesher::param {

// The two sides of a chain of edges on a triangle mesh.
//  - every chain edge gets its two corners as a pair (a, b); for the first
//    edge in corner order (a border edge has no second corner, an edge with
//    more than two corners none at all),
//  - side A is handed on along the chain: from corner a of the edge before,
//    the triangles around the chain point between the two edges are walked
//    until the next chain edge is reached; its corner in the triangle the
//    walk arrives in becomes a, the other one b,
//  - the walk fails at a missing edge and at any edge without exactly two
//    corners, so a border edge can only be the first edge of a chain; a walk
//    that comes back to where it began fails as well.
class ChainSides {
public:
    // Per chain edge its corner on side A and on side B (-1: none).
    struct CornerPairs {
        std::vector<int> a;
        std::vector<int> b;
    };
    // `points`: the chain points, one more than `edges`. Empty when a walk failed.
    [[nodiscard]] static std::optional<CornerPairs> corners(const mesh::Adjacency& adjacency,
                                                            std::span<const int> edges, std::span<const int> points);

    // Sides of one chain on the cut copy.
    struct Sides {
        bool found = false;
        // Quarter turns left across the first chain edge, counted from side A.
        int jump = 0;
        // Per chain point the copy point on side A and on side B: the copy
        // point of the corner of the pair of the edge that starts there (for
        // the last point the last edge) where that corner lies at the chain
        // point, else of its next corner. A point whose pair lacks a corner
        // is left out of both lists.
        std::vector<int> a;
        std::vector<int> b;
    };
    // `adjacency`: the neighbourhood of the mesh the chain lies on, `copy`
    // its cut copy, `residualJump`: per edge of the mesh the quarter turns
    // left across it. `edges` holds at least one edge. Without corner pairs
    // nothing is found.
    [[nodiscard]] static Sides of(const mesh::Adjacency& adjacency, const field::CutCopy& copy,
                                  std::span<const int> residualJump, std::span<const int> edges,
                                  std::span<const int> points);
};

} // namespace areko::remesher::param
