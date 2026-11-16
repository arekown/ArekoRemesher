#pragma once

#include "adjacency.h"
#include "bit_array.h"
#include "csr.h"

#include <span>
#include <vector>

namespace areko::remesher::mesh {

// Splits a set of edges into point chains.
//  - degree d(v) = number of chosen edges at v; points are visited in the
//    order they first appear as end points of the chosen edges (ascending),
//  - chains start at ends (d = 1, last one first), then at branch points
//    (d >= 3, d-1 starts each, taken from the back), then once more at the
//    branch points after refilling the start stack (its first entry stays),
//    finally (or at once without branch points) at the first end point of the
//    smallest free edge,
//  - a walk follows the first free edge in the point's edge order and stops at
//    the first point with d != 2; a closed loop repeats its start point.
// Walks without an edge are dropped. Along with the points of every chain its
// edges are found (row r: edges between consecutive points of chain r).
// Separators: each listed point with d = 2 (in list order, after the fallback
// start of a pure loop) becomes a branch point with d = 3, so chains end and
// start there; other listed points change nothing.
class EdgeChains {
public:
    // Where a caller departs from the rules above.
    struct Options {
        std::span<const int> separators;
        // The border points of the mesh separate like listed separators, after them.
        bool borderPointsSeparate = false;
        // Starts follow the chains: the starts at branch points are not laid
        // out beforehand but added as chains end there (d-1 each time a chain
        // ends at a point with d >= 3) and taken in the order they were added;
        // when they run out, the starts of all branch points are added once.
        bool startsFollowChains = false;
    };

    struct Result {
        Csr points; // chain -> its points
        Csr edges;  // chain -> the edges between its consecutive points
        // The nodes of the chains: the ends in the order the points were
        // visited, then the branch points (visiting order; for pure loops the
        // first visited point, which separates nothing), then the separators
        // and border points that became branch points. The start of pure
        // loops is listed once, also when it became a branch point.
        std::vector<int> nodes;
    };

    // The points of the chains; with `chainEdgeRows` their edges as well.
    [[nodiscard]] static Csr build(const Adjacency& adjacency, const BitArray& chosen, Csr* chainEdgeRows = nullptr,
                                   std::span<const int> separators = {});
    [[nodiscard]] static Result build(const Adjacency& adjacency, const BitArray& chosen, const Options& options);
    // The same for a list of edges: the points are visited in the order of
    // the list instead of ascending edges; an edge listed again counts once.
    [[nodiscard]] static Result build(const Adjacency& adjacency, std::span<const int> edgeList,
                                      const Options& options);

private:
    // `order`: the chosen edges in the order they are taken up, each once, or
    // empty for ascending.
    [[nodiscard]] static Result run(const Adjacency& adjacency, const BitArray& chosen, std::span<const int> order,
                                    const Options& options);
};

} // namespace areko::remesher::mesh
