#pragma once

#include "cut_graph.h"

#include "mesh/adjacency.h"

#include <cstddef>
#include <optional>

namespace areko::remesher::field {

// From which neighbour a triangle takes over values when values per triangle
// and combed direction are smoothed, and which of them:
//  - the neighbour is the triangle across the edge of a corner; a corner
//    without edge has none, and neither has an edge that does not join
//    exactly two corners,
//  - across an edge with an odd number of quarter turns left between the
//    combed crosses the two directions of the neighbour change places,
//  - across a curve edge that limits, only the direction along the curve
//    passes: the first one when the curve follows it at that corner, else the
//    second.
// The sums stay with the caller.
class SmoothingNeighbours {
public:
    // Which curve edges limit what passes.
    enum class Limit {
        EveryCurve, // every curve edge
        Creases,    // creases only (also border edges and edges on the mirror plane)
    };
    // Which of the two directions of a neighbour count.
    struct Take {
        bool u = true;
        bool v = true;
    };
    struct Neighbour {
        int triangle = -1;
        bool odd = false; // its two directions change places
        Take take;
    };

    SmoothingNeighbours(const mesh::Adjacency& adjacency, const CutGraph& cuts, Limit limit) noexcept
        : adjacency_(adjacency), cuts_(cuts), limit_(limit)
    {
    }

    // What passes edge `edge`, seen from its corner `corner`.
    [[nodiscard]] Take take(int corner, int edge) const
    {
        const std::uint8_t kind = cuts_.curveEdgeKind[static_cast<std::size_t>(edge)];
        if (limit_ == Limit::Creases ? kind != curve_edge::crease : kind == curve_edge::none) {
            return {};
        }
        const bool alongU = cuts_.curveFollowsU.testBounded(static_cast<std::size_t>(corner));
        return {alongU, !alongU};
    }

    // The neighbour across the edge of `corner`, if there is one.
    [[nodiscard]] std::optional<Neighbour> across(int corner) const
    {
        const int e = adjacency_.cornerEdge(corner);
        if (e < 0) {
            return std::nullopt;
        }
        const int other = adjacency_.otherCorner(e, corner);
        if (other < 0) {
            return std::nullopt;
        }
        const bool odd = (cuts_.residualJump[static_cast<std::size_t>(e)] & 1) != 0;
        return Neighbour{other / 3, odd, take(corner, e)};
    }

private:
    const mesh::Adjacency& adjacency_;
    const CutGraph& cuts_;
    Limit limit_;
};

} // namespace areko::remesher::field
