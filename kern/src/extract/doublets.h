#pragma once

#include "mesh/mesh.h"

#include <span>

namespace areko::remesher::extract {

// Dissolves doublets of the raw mesh: up to three rounds; each round lists
// both edges of every point with exactly two edges and exactly two corners
// (points ascending) and dissolves them with nearly straight points
// (EdgeDissolve). A round without such a point ends the work. Border points
// are not treated apart.
class Doublets {
public:
    static constexpr int rounds = 3;

    // Returns the number of rounds that dissolved something.
    static int run(mesh::Mesh& mesh);

    // Calls `visit(edges, corners)` for every point with exactly two edges
    // and exactly two corners, points ascending. Which of them count as
    // doublets is the rule of the caller.
    template <class Visit>
    static void forEachPoint(const mesh::Mesh& mesh, Visit&& visit)
    {
        const mesh::Adjacency& adjacency = mesh.adjacency();
        for (int v = 0; v < mesh.pointCount(); ++v) {
            const std::span<const int> edges = adjacency.pointEdges(v);
            const std::span<const int> corners = adjacency.pointCorners(v);
            if (edges.size() == 2 && corners.size() == 2) {
                visit(edges, corners);
            }
        }
    }
};

} // namespace areko::remesher::extract
