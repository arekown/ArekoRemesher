#include "doublets.h"

#include "edge_dissolve.h"

#include <vector>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

int Doublets::run(Mesh& mesh)
{
    int done = 0;
    std::vector<int> edges;
    for (int round = 0; round < rounds; ++round) {
        edges.clear();
        forEachPoint(mesh, [&](std::span<const int> pointEdges, std::span<const int>) {
            edges.insert(edges.end(), pointEdges.begin(), pointEdges.end());
        });
        if (edges.empty()) {
            break;
        }
        EdgeDissolve::run(mesh, edges, true);
        ++done;
    }
    return done;
}

} // namespace areko::remesher::extract
