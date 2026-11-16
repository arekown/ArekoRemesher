#include "dual_spanning_tree.h"

#include "centre_path_search.h"

#include <vector>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

BitArray DualSpanningTree::freeEdges(const Mesh& mesh, const BitArray& anchors, const BitArray& mirrorHalf)
{
    const Adjacency& adjacency = mesh.adjacency();
    const std::size_t faceCount = static_cast<std::size_t>(mesh.faceCount());
    const std::size_t edgeCount = static_cast<std::size_t>(adjacency.edgeCount());
    if (faceCount == 0) {
        return BitArray(edgeCount);
    }
    BitArray free(edgeCount, true);
    std::vector<int> parent(faceCount, -1);
    CentrePathSearch search(mesh);
    anchors.forEachSet([&](std::size_t t) {
        if (!mirrorHalf.testBounded(t)) {
            search.start(static_cast<int>(t));
        }
    });
    if (search.empty()) {
        search.start(0);
    }
    // The parent of a triangle is the edge its shortest path arrives across.
    search.run([&](int n, int, int e) { parent[static_cast<std::size_t>(n)] = e; });
    for (const int e : parent) {
        if (e != -1) {
            free.set(static_cast<std::size_t>(e), false);
        }
    }
    for (std::size_t e = 0; e < edgeCount; ++e) {
        const auto faces = adjacency.innerEdgeFaces(static_cast<int>(e));
        if (!faces || (anchors.testBounded(static_cast<std::size_t>((*faces)[0])) &&
                       anchors.testBounded(static_cast<std::size_t>((*faces)[1])))) {
            free.set(e, false);
        }
    }
    return free;
}

} // namespace areko::remesher::field
