#pragma once

#include "prep_state.h"

#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Cleans and repairs the working mesh of the preparation:
//  - mode 0 first removes repeated faces,
//  - joins the end points of very short edges (limits from the bounding box
//    diagonal d: (1e-8 d)^2, at least 9.999999e-19, and (5e-7 d)^2),
//  - repairs the topology,
//  - closes border loops of up to four points whose perimeter is below
//    0.25 L' or whose area is below L'^2 / 50 (L' from the target edge length
//    or the area with at least 3000 quads), then repairs again,
//  - dissolves points between two edges and drops unused points,
//  - triangulates again when the mesh is no triangle mesh any more.
// Writes the repair count of the state; L' is always derived here.
class MeshCleaner {
public:
    static void run(mesh::Mesh& mesh, int mode, PrepState& state);
};

} // namespace areko::remesher::prep
