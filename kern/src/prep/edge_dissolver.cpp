#include "edge_dissolver.h"

#include "mesh/face_rebuild.h"

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

int EdgeDissolver::run(Mesh& mesh)
{
    const Adjacency& adjacency = mesh.adjacency();
    std::vector<char> dissolve(static_cast<std::size_t>(mesh.pointCount()), 0);
    int marked = 0;
    for (int v = 0; v < mesh.pointCount(); ++v) {
        if (adjacency.pointEdges(v).size() == 2 && adjacency.pointCorners(v).size() == 2) {
            dissolve[static_cast<std::size_t>(v)] = 1;
            ++marked;
        }
    }
    const Topology& topology = mesh.topology();
    if (marked == 0 || topology.mode() != FaceMode::Polygons) {
        return 0;
    }
    FaceRebuild rebuild(mesh, topology.cornerVertex());
    std::vector<int> kept;
    for (int f = 0; f < topology.faceCount(); ++f) {
        kept.clear();
        for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
            if (!dissolve[static_cast<std::size_t>(topology.vertex(c))]) {
                kept.push_back(c);
            }
        }
        rebuild.keep(f, kept);
    }
    rebuild.apply(mesh, FaceMode::Polygons);
    return marked;
}

} // namespace areko::remesher::prep
