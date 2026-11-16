#include "normal_split_edges.h"

#include "mesh/sharp_edges.h"

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

void NormalSplitEdges::apply(Mesh& mesh)
{
    const auto& normals = mesh.customNormals();
    if (!normals || normals->cornerNormal.size() != static_cast<std::size_t>(mesh.cornerCount())) {
        return;
    }
    const auto& index = normals->cornerNormal;
    const Adjacency& adjacency = mesh.adjacency();
    BitArray sharp = SharpEdges::extendableEdgeBits(mesh);
    for (int e = 0; e < adjacency.edgeCount(); ++e) {
        const auto corners = adjacency.edgeCorners(e);
        if (corners.size() != 2) {
            continue;
        }
        const int c0 = corners[0];
        const int c1 = corners[1];
        const int n0 = adjacency.next(c0);
        const int n1 = adjacency.next(c1);
        if (adjacency.vertex(c0) == adjacency.vertex(n1) &&
            (index[static_cast<std::size_t>(c0)] != index[static_cast<std::size_t>(n1)] ||
             index[static_cast<std::size_t>(n0)] != index[static_cast<std::size_t>(c1)])) {
            sharp.set(static_cast<std::size_t>(e));
        }
    }
    SharpEdges::setEdgeBits(mesh, std::move(sharp));
}

} // namespace areko::remesher::prep
