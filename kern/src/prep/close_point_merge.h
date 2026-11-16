#pragma once

#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Joins the end points of very short edges. An edge is a candidate when its
// squared length (float) is at most `nearSquared`, or, on a triangle mesh with
// `smallSquared` > 0, at most `smallSquared` while none of its triangles has
// a doubled area above 1e-8 times the total area (at least 1e-12; a NaN area
// counts as small). Without both rules only edges with exactly equal end
// points count. Connected candidates form groups; every group joins onto one
// existing point (its position stays): its first border point in the order
// the group lists its points, else the end point [0] of its first edge. Faces
// at moved points lose corners that repeat the next corner and vanish below
// three corners. Points stay (the moved ones become unused); point maps,
// face maps and sharp edges follow.
class ClosePointMerge {
public:
    // Returns the number of points moved onto another one.
    static int apply(mesh::Mesh& mesh, float nearSquared, float smallSquared);
};

} // namespace areko::remesher::prep
