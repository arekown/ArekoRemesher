#pragma once

#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Edges where the given corner normals split become sharp edges: an inner edge
// with two corners of opposite orientation is a split when the normal index
// differs at one of its two end points. The splits are added to the sharp
// edges already present; without given normals nothing changes.
class NormalSplitEdges {
public:
    static void apply(mesh::Mesh& mesh);
};

} // namespace areko::remesher::prep
