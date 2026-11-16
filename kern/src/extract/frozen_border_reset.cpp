#include "frozen_border_reset.h"

#include "mesh/boundary_loops.h"

namespace areko::remesher::extract {

void FrozenBorderReset::run(mesh::Mesh& output, const mesh::Mesh& triangles,
                            const std::vector<mesh::FreezeMark>& marks, const std::vector<int>& originalToOutput)
{
    const mesh::BoundaryLoops loops(triangles);
    std::vector<int> moved;
    std::vector<mesh::Vec3f> positions;
    for (const int v : loops.points().entries) {
        const auto i = static_cast<std::size_t>(v);
        const int out = originalToOutput[i];
        if (marks[i].value == 1 && out != -1) {
            moved.push_back(out);
            positions.push_back(triangles.points()[i]);
        }
    }
    if (!moved.empty()) {
        output.movePoints(moved, positions);
    }
}

} // namespace areko::remesher::extract
