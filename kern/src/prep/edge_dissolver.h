#pragma once

#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Dissolves points that only join two edges: a point with exactly two edges
// and exactly two corners leaves its faces, so its two edges become one.
// Faces with fewer than three corners left vanish. Works on polygon meshes
// only; points stay (the dissolved ones become unused).
class EdgeDissolver {
public:
    // Returns the number of dissolved points (0 without polygon mode).
    static int run(mesh::Mesh& mesh);
};

} // namespace areko::remesher::prep
