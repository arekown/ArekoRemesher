#include "link_path.h"

#include <limits>

namespace areko::remesher::mesh {

LinkPath shortestLink(const Adjacency& adjacency, std::span<const Vec3f> points, int a, int b)
{
    LinkPath path;
    const int direct = adjacency.edgeBetween(a, b);
    if (direct >= 0) {
        path.first = direct;
        return path;
    }
    const auto edgesB = adjacency.pointEdges(b);
    const auto neighboursB = adjacency.pointNeighbours(b);
    float best = std::numeric_limits<float>::max();
    for (std::size_t i = 0; i < edgesB.size(); ++i) {
        const int c = neighboursB[i];
        const int viaEdge = adjacency.edgeBetween(a, c);
        if (viaEdge < 0) {
            continue;
        }
        const Vec3f pc = points[static_cast<std::size_t>(c)];
        const float sum = distance(pc, points[static_cast<std::size_t>(a)]) + distance(pc, points[static_cast<std::size_t>(b)]);
        if (sum < best) {
            best = sum;
            path.first = viaEdge;
            path.second = edgesB[i];
            path.via = c;
        }
    }
    return path;
}

} // namespace areko::remesher::mesh
