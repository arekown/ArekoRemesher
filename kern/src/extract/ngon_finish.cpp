#include "ngon_finish.h"

#include "doublets.h"

#include "mesh/bit_array.h"

#include <algorithm>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

std::vector<int> DoubletEdges::find(const Mesh& mesh)
{
    const Adjacency& adjacency = mesh.adjacency();
    const Topology& topology = mesh.topology();
    BitArray chosen(static_cast<std::size_t>(adjacency.edgeCount()));
    if (topology.mode() != FaceMode::Triangles) {
        Doublets::forEachPoint(mesh, [&](std::span<const int> edges, std::span<const int> corners) {
            if (adjacency.isBoundaryEdge(edges[0]) || adjacency.isBoundaryEdge(edges[1])) {
                return;
            }
            if (topology.mode() == FaceMode::Polygons &&
                (topology.faceSize(adjacency.cornerFace(corners[0])) != 4 ||
                 topology.faceSize(adjacency.cornerFace(corners[1])) != 4)) {
                return;
            }
            chosen.set(static_cast<std::size_t>(edges.front()));
            chosen.set(static_cast<std::size_t>(edges.back()));
        });
    }
    std::vector<int> out;
    chosen.forEachSet([&](std::size_t e) { out.push_back(static_cast<int>(e)); });
    return out;
}

void NewPointSmoother::run(Mesh& mesh, int first, int iterations, float lambda)
{
    const int count = mesh.pointCount();
    if (count == first) {
        return;
    }
    const Adjacency& adjacency = mesh.adjacency();
    std::vector<Vec3f> read = mesh.points();
    // Start from the current place, so a point without edges keeps it.
    std::vector<Vec3f> moved(read.begin() + first, read.end());
    for (int iteration = 0; iteration < iterations; ++iteration) {
        for (int v = first; v < count; ++v) {
            const auto neighbours = adjacency.pointNeighbours(v);
            if (neighbours.empty()) {
                continue;
            }
            Vec3f sum{};
            for (const int w : neighbours) {
                sum = sum + read[static_cast<std::size_t>(w)];
            }
            const float share = lambda / static_cast<float>(neighbours.size());
            moved[static_cast<std::size_t>(v - first)] = read[static_cast<std::size_t>(v)] * (1.0f - lambda) + sum * share;
        }
        std::copy(moved.begin(), moved.end(), read.begin() + first);
    }
    std::vector<int> indices(moved.size());
    for (std::size_t i = 0; i < indices.size(); ++i) {
        indices[i] = first + static_cast<int>(i);
    }
    mesh.movePoints(indices, moved);
}

} // namespace areko::remesher::extract
