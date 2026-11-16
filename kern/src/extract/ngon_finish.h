#pragma once

#include "mesh/mesh.h"

#include <vector>

namespace areko::remesher::extract {

// Edges of the doublets of a mesh: points with exactly two edges and two
// corners, neither edge on the border; in polygon mode both faces must be
// quads, in quad mode every such point counts, in triangle mode none.
// Both edges of each such point, ascending, each once.
class DoubletEdges {
public:
    [[nodiscard]] static std::vector<int> find(const mesh::Mesh& mesh);
};

// Moves the points from `first` on towards the average of their edge
// neighbours: per iteration (all read from the previous one)
// p' = (1 - lambda) p + sum * (lambda / count). A point without edges stays.
class NewPointSmoother {
public:
    static void run(mesh::Mesh& mesh, int first, int iterations, float lambda);
};

} // namespace areko::remesher::extract
