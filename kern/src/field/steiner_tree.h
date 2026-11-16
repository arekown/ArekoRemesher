#pragma once

#include "mesh/adjacency.h"
#include "mesh/bit_array.h"

#include <span>
#include <vector>

namespace areko::remesher::field {

// Tree of mesh edges that joins a set of sources along cheap paths.
//  - every point source and every line source (a set of points that already
//    hang together, e.g. a border loop) starts a region; line sources that
//    touch an earlier source are joined with it beforehand,
//  - the regions grow over the edges by their weights (first in, first out
//    with later improvements; border edges are not walked, border points not
//    entered),
//  - for every pair of regions that meet, the cheapest edge between them is a
//    bridge: the distances of its two ends plus its weight,
//  - the bridges are taken cheapest first as long as they join two regions
//    that are not joined yet; each taken bridge contributes its edge and the
//    paths from both ends back to their sources.
class SteinerTree {
public:
    struct Sources {
        std::vector<int> points;             // a point named twice counts once
        std::vector<std::vector<int>> lines; // empty lines are skipped
    };

    // One bit per edge: the edges of the tree.
    [[nodiscard]] static mesh::BitArray edges(const mesh::Adjacency& adjacency, std::span<const float> weights,
                                              const Sources& sources);
};

} // namespace areko::remesher::field
