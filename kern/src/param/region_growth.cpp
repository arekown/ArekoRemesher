#include "region_growth.h"

#include <cfloat>
#include <span>

namespace areko::remesher::param {

RegionGrowth::RegionGrowth(const SystemContext& context, const UvOutput::Coordinates& uv, int coordinate)
    : copy_(context.copy().adjacency()),
      coordinate_(coordinate),
      chainRegions_(context.chains().curves.count()),
      distance_(uv.size(), FLT_MAX),
      region_(uv.size(), -1),
      root_(uv.size(), -1),
      previousEdge_(uv.size(), -1),
      chainEdges_(static_cast<std::size_t>(copy_.edgeCount()))
{
    util::IndexQueue queue(uv.size());
    addChains(context, queue);
    addSingularPoints(context, queue);
    grow(uv, queue);
}

void RegionGrowth::makeSource(int w, int region, util::IndexQueue& queue)
{
    const auto i = static_cast<std::size_t>(w);
    distance_[i] = 0.0f;
    region_[i] = region;
    root_[i] = w;
    previousEdge_[i] = -1;
    queue.pushBack(w);
}

// The curve chains along which the coordinate stays constant: the first grid
// coordinate along a chain that follows the second combed direction and the
// other way round. The first corner of the first edge of a chain decides.
void RegionGrowth::addChains(const SystemContext& context, util::IndexQueue& queue)
{
    const mesh::Csr& chainEdges = context.chains().curves.edges;
    const mesh::BitArray& followsFirst = context.graph().curveFollowsU;
    for (int k = 0; k < chainRegions_; ++k) {
        const std::span<const int> edges = chainEdges.row(k);
        if (edges.empty()) {
            continue;
        }
        const int corner = copy_.edgeCorners(edges.front())[0];
        const bool alongFirst = followsFirst.testBounded(static_cast<std::size_t>(corner));
        if (alongFirst != (coordinate_ == 1)) {
            continue;
        }
        for (const int e : edges) {
            chainEdges_.set(static_cast<std::size_t>(e));
            for (const int w : copy_.edgeVertices(e)) {
                if (!context.mirrored(w)) {
                    makeSource(w, k, queue);
                }
            }
        }
    }
}

// Every copy of a singular point that no chain has taken.
void RegionGrowth::addSingularPoints(const SystemContext& context, util::IndexQueue& queue)
{
    const std::vector<int>& points = context.chains().singularities.points;
    for (std::size_t j = 0; j < points.size(); ++j) {
        const int region = chainRegions_ + coordinate_ + 2 * static_cast<int>(j);
        for (const int w : context.copies(points[j])) {
            if (region_[static_cast<std::size_t>(w)] < 0 && !context.mirrored(w)) {
                makeSource(w, region, queue);
            }
        }
    }
}

// Shortest paths from all sources at once: a point whose distance got
// smaller is looked at again, first in, first out. Only a strictly smaller
// distance replaces the one a point has.
void RegionGrowth::grow(const UvOutput::Coordinates& uv, util::IndexQueue& queue)
{
    while (!queue.empty()) {
        const int v = queue.pop();
        const auto from = static_cast<std::size_t>(v);
        for (const int e : copy_.pointEdges(v)) {
            const int o = copy_.otherEnd(e, v);
            const auto to = static_cast<std::size_t>(o);
            const float reached = gridDistance(uv, v, o) + distance_[from];
            if (reached < distance_[to]) {
                distance_[to] = reached;
                region_[to] = region_[from];
                root_[to] = root_[from];
                previousEdge_[to] = e;
                queue.pushBack(o);
            }
        }
    }
}

std::vector<int> RegionGrowth::pathToRoot(int w) const
{
    std::vector<int> path{w};
    for (int e = previousEdge_[static_cast<std::size_t>(w)]; e >= 0; e = previousEdge_[static_cast<std::size_t>(w)]) {
        w = copy_.otherEnd(e, w);
        path.push_back(w);
    }
    return path;
}

} // namespace areko::remesher::param
