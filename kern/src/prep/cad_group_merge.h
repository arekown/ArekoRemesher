#pragma once

#include "group_union.h"
#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Merges narrow strips of CAD faces into their neighbours. In up to five
// passes the borders between face groups (and the mesh boundary) are split
// into chains; a short chain whose both ends meet two other chains at about
// right angles that run on in a straight line is a strip edge. All strip edges
// of a pass are removed from the borders at once and the face groups become
// the new components. Another pass follows while a long chain touches a chain
// with an end at a removed strip.
class CadGroupMerge {
public:
    // The restore table follows the renumbering through the old group numbers.
    static void apply(mesh::Mesh& mesh, GroupRestoreTable& table);
};

} // namespace areko::remesher::prep
