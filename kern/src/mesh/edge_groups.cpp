#include "edge_groups.h"

#include <vector>

namespace areko::remesher::mesh {

Csr EdgeGroups::byConnectivity(const Adjacency& adjacency, std::span<const int> edges)
{
    std::vector<char> listed(static_cast<std::size_t>(adjacency.edgeCount()), 0);
    for (const int e : edges) {
        listed[static_cast<std::size_t>(e)] = 1;
    }
    std::vector<char> taken(listed.size(), 0);
    std::vector<char> onStack(static_cast<std::size_t>(adjacency.pointCount()), 0);
    Csr groups;
    std::vector<int> stack;
    for (const int seed : edges) {
        if (taken[static_cast<std::size_t>(seed)]) {
            continue;
        }
        stack.assign(1, adjacency.edgeVertices(seed)[0]);
        onStack[static_cast<std::size_t>(stack.back())] = 1;
        while (!stack.empty()) {
            const int v = stack.back();
            stack.pop_back();
            onStack[static_cast<std::size_t>(v)] = 0;
            for (const int e : adjacency.pointEdges(v)) {
                if (!listed[static_cast<std::size_t>(e)] || taken[static_cast<std::size_t>(e)]) {
                    continue;
                }
                taken[static_cast<std::size_t>(e)] = 1;
                groups.entries.push_back(e);
                const int w = adjacency.otherEnd(e, v);
                if (!onStack[static_cast<std::size_t>(w)]) {
                    onStack[static_cast<std::size_t>(w)] = 1;
                    stack.push_back(w);
                }
            }
        }
        groups.starts.push_back(static_cast<int>(groups.entries.size()));
    }
    return groups;
}

} // namespace areko::remesher::mesh
