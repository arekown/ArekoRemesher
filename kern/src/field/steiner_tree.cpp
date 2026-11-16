#include "steiner_tree.h"

#include "util/disjoint_sets.h"
#include "util/introsort.h"
#include "util/reach_queue.h"

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <unordered_map>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

BitArray SteinerTree::edges(const Adjacency& adjacency, std::span<const float> weights, const Sources& sources)
{
    const std::size_t pointCount = static_cast<std::size_t>(adjacency.pointCount());
    const int edgeCount = adjacency.edgeCount();
    BitArray tree(static_cast<std::size_t>(edgeCount));
    util::ReachQueue search(pointCount);
    std::vector<int> region(pointCount, -1);
    std::vector<int> cameBy(pointCount, -1);
    int regionCount = 0;
    auto claim = [&](int v) {
        region[static_cast<std::size_t>(v)] = regionCount;
        search.start(v);
    };
    auto isFree = [&](int v) { return !search.reached(v); };

    for (const int v : sources.points) {
        if (isFree(v)) {
            claim(v);
            ++regionCount;
        }
    }
    // A line that touches earlier sources forms a group with their regions.
    std::vector<std::vector<int>> groups;
    for (const std::vector<int>& line : sources.lines) {
        if (line.empty()) {
            continue;
        }
        const bool anyFree = std::any_of(line.begin(), line.end(), isFree);
        std::vector<int> group;
        for (const int v : line) {
            if (!isFree(v)) {
                group.push_back(region[static_cast<std::size_t>(v)]);
            }
        }
        if (!group.empty()) {
            if (anyFree) {
                group.push_back(regionCount);
            }
            groups.push_back(std::move(group));
        }
        for (const int v : line) {
            if (isFree(v)) {
                claim(v);
            }
        }
        if (anyFree) {
            ++regionCount;
        }
    }

    // Shortest paths from all sources at once: first in, first out, a point
    // improved later is queued again.
    BitArray borderPoint(pointCount);
    for (const int v : adjacency.boundaryPoints()) {
        borderPoint.set(static_cast<std::size_t>(v));
    }
    search.run([&](int v) {
        for (const int e : adjacency.pointEdges(v)) {
            if (adjacency.isBoundaryEdge(e)) {
                continue;
            }
            const int u = adjacency.otherEnd(e, v);
            if (borderPoint.test(static_cast<std::size_t>(u))) {
                continue;
            }
            if (search.offer(u, search.reach(v) + weights[static_cast<std::size_t>(e)])) {
                region[static_cast<std::size_t>(u)] = region[static_cast<std::size_t>(v)];
                cameBy[static_cast<std::size_t>(u)] = e;
            }
        }
    });

    // Bridges: per pair of regions the cheapest edge between them, in the
    // order the pairs first appear.
    struct Bridge {
        int low;
        int high;
        int edge;
        float cost;
    };
    std::vector<Bridge> bridges;
    std::unordered_map<std::uint64_t, std::size_t> bridgeOfPair;
    for (int e = 0; e < edgeCount; ++e) {
        if (adjacency.isBoundaryEdge(e)) {
            continue;
        }
        const auto ends = adjacency.edgeVertices(e);
        const int a = region[static_cast<std::size_t>(ends[0])];
        const int b = region[static_cast<std::size_t>(ends[1])];
        if (a == -1 || b == -1 || a == b) {
            continue;
        }
        const float cost = (search.reach(ends[0]) + search.reach(ends[1])) + weights[static_cast<std::size_t>(e)];
        const int low = std::min(a, b);
        const int high = std::max(a, b);
        const std::uint64_t key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(low)) << 32) |
                                  static_cast<std::uint32_t>(high);
        const auto [found, isNew] = bridgeOfPair.try_emplace(key, bridges.size());
        if (isNew) {
            bridges.push_back({low, high, e, cost});
        } else if (cost < bridges[found->second].cost) {
            bridges[found->second].edge = e;
            bridges[found->second].cost = cost;
        }
    }
    std::vector<int> order(bridges.size());
    std::iota(order.begin(), order.end(), 0);
    util::Introsort::sort(order.begin(), order.end(), [&](int a, int b) {
        return bridges[static_cast<std::size_t>(a)].cost < bridges[static_cast<std::size_t>(b)].cost;
    });

    // Which regions are joined; the regions of a group touch from the start.
    util::DisjointSets regions(static_cast<std::size_t>(regionCount));
    for (const std::vector<int>& group : groups) {
        regions.joinAll(group);
    }
    // Cheapest bridges first, each only between regions not joined yet; a
    // taken bridge marks its edge and both paths back to the sources.
    auto markPath = [&](int v) {
        for (int e = cameBy[static_cast<std::size_t>(v)]; e != -1; e = cameBy[static_cast<std::size_t>(v)]) {
            tree.set(static_cast<std::size_t>(e));
            v = adjacency.otherEnd(e, v);
        }
    };
    for (const int i : order) {
        const Bridge& bridge = bridges[static_cast<std::size_t>(i)];
        if (regions.joined(bridge.low, bridge.high)) {
            continue;
        }
        regions.join(bridge.low, bridge.high);
        tree.set(static_cast<std::size_t>(bridge.edge));
        const auto ends = adjacency.edgeVertices(bridge.edge);
        markPath(ends[0]);
        markPath(ends[1]);
    }
    return tree;
}

} // namespace areko::remesher::field
