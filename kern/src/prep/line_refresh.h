#pragma once

#include "mesh/mesh.h"
#include "prep_state.h"

namespace areko::remesher::prep {

// Renews the lines from the sharp edges and improves flat, acute triangles of
// the triangle mesh: up to 5 runs of edge flips, then the treatment of tiny
// triangles. Face group borders, the mesh border and the line edges are
// protected; the lines follow the flips.
class LineRefresh {
public:
    static void run(mesh::Mesh& mesh, PrepState& state);
};

} // namespace areko::remesher::prep
