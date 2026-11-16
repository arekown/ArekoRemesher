#pragma once

#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Moves tiny needle and cap triangles into the face group of a neighbour, so
// that group borders do not run through them; sharp edge bits between the
// triangle and its new group are cleared. Positions and faces stay unchanged.
//  - "bad" (decided once): doubled area below 0.2 L^2 and one corner angle
//    below 0.3 degrees or two below 2 degrees,
//  - up to 10 rounds over all bad triangles in ascending order, a further
//    round only after a change,
//  - needs a triangle mesh with face groups.
class DegenerateGroups {
public:
    static void run(mesh::Mesh& mesh, float length);
};

} // namespace areko::remesher::prep
