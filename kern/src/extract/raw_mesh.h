#pragma once

#include "output_points.h"

#include "mesh/mesh.h"

#include <vector>

namespace areko::remesher::extract {

// The raw quad mesh of one part: the polygon mesh built from the closed faces
// of the trace (object coordinates, unused points included), what every
// output point comes from, and per point of the triangle mesh its output
// point (-1 for none).
struct RawMesh {
    mesh::Mesh mesh;
    OutputPoints points;
    std::vector<int> originalToOutput;
    // Edges that follow a crease curve (one bit per edge; set when the
    // features are transferred, empty before and without curves).
    mesh::BitArray curveEdges;
};

} // namespace areko::remesher::extract
