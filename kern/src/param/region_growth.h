#pragma once

#include "mesh/bit_array.h"
#include "param/system_context.h"
#include "param/uv_output.h"
#include "util/index_queue.h"

#include <cmath>
#include <cstddef>
#include <vector>

namespace areko::remesher::param {

// Distance of two copy points in the grid coordinates (float).
[[nodiscard]] inline float gridDistance(const UvOutput::Coordinates& uv, int a, int b)
{
    const float du = uv[static_cast<std::size_t>(a)][0] - uv[static_cast<std::size_t>(b)][0];
    const float dv = uv[static_cast<std::size_t>(a)][1] - uv[static_cast<std::size_t>(b)][1];
    return std::sqrt(du * du + dv * dv);
}

// Regions of one grid coordinate on the cut copy of a part. The sources are
// the curve chains along which this coordinate stays constant and the
// singular points; every point of the copy belongs to the source it is
// nearest to, measured along the edges in the grid coordinates of the latest
// solution. Between the sources of two neighbouring regions a whole number
// of grid lines has to fit.
//  - region of curve chain k: k; region of singular point j: number of curve
//    chains + coordinate + 2j,
//  - a point of the mirrored half is never a source but is reached,
//  - of two chains through a point the later one has it,
//  - with equal distances the source that got there first keeps the point.
class RegionGrowth {
public:
    // `coordinate`: 0 for the first grid coordinate, 1 for the second. `uv`
    // holds the coordinates of every copy point.
    RegionGrowth(const SystemContext& context, const UvOutput::Coordinates& uv, int coordinate);

    [[nodiscard]] int coordinate() const noexcept { return coordinate_; }
    // Regions below this number belong to curve chains.
    [[nodiscard]] int chainRegions() const noexcept { return chainRegions_; }

    // Region of copy point w (-1: not reached), the source point its path
    // starts at and the length of that path.
    [[nodiscard]] int region(int w) const { return region_[static_cast<std::size_t>(w)]; }
    [[nodiscard]] int root(int w) const { return root_[static_cast<std::size_t>(w)]; }
    [[nodiscard]] float distance(int w) const { return distance_[static_cast<std::size_t>(w)]; }
    // Whether copy edge e lies on a chain that is a source.
    [[nodiscard]] bool onChain(int e) const { return chainEdges_.testBounded(static_cast<std::size_t>(e)); }

    // The points from w back to the point its path starts at, w first.
    [[nodiscard]] std::vector<int> pathToRoot(int w) const;

private:
    // `queue` holds the points whose neighbours are still to be looked at.
    void makeSource(int w, int region, util::IndexQueue& queue);
    void addChains(const SystemContext& context, util::IndexQueue& queue);
    void addSingularPoints(const SystemContext& context, util::IndexQueue& queue);
    void grow(const UvOutput::Coordinates& uv, util::IndexQueue& queue);

    const mesh::Adjacency& copy_;
    int coordinate_;
    int chainRegions_;
    std::vector<float> distance_;
    std::vector<int> region_;
    std::vector<int> root_;
    std::vector<int> previousEdge_; // edge towards the source, -1 at a source
    mesh::BitArray chainEdges_;
};

} // namespace areko::remesher::param
