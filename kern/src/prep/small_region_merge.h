#pragma once

#include "prep_state.h"

#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Merges small face group regions into a large neighbour (triangle meshes
// with face groups; only groups and sharp edge bits change):
//  - regions are the parts between face group borders and sharp edges, with
//    their doubled area; the threshold is 2 A / (15 N) for the clamped quad
//    count N (0: 5000, 1000 .. 10000),
//  - one pass over the regions in order; areas, regions and borders stay as
//    they were at the start, the groups change at once,
//  - a small region chooses a large neighbour group (largest shared border
//    for up to three neighbours and few faces, otherwise the only large one)
//    and merges when its shape is thin or jagged; the border measure carries
//    over from region to region (half of previous plus own sum, from 0),
//  - merging clears the sharp bits towards faces of the new group.
class SmallRegionMerge {
public:
    static void run(mesh::Mesh& mesh, const PrepState& state);
};

} // namespace areko::remesher::prep
