#pragma once

#include "mesh/mesh.h"
#include "prep/group_union.h"

namespace areko::remesher::merge {

// Gives the result the materials, face groups and smoothing groups of the
// host back. The face groups of the result are region numbers of the
// preparation; nothing happens when the result has no face groups or fewer
// values than faces. Per face with region g:
//  - material: the material of region g, or g itself when the face groups
//    were a copy of the materials, otherwise unchanged,
//  - face group: the host face group of region g,
//  - smoothing groups: those of region g,
// each only when the table has such values; a region beyond the table
// gives 0.
class GroupRestore {
public:
    static void run(mesh::Mesh& mesh, const prep::GroupRestoreTable& table);
};

} // namespace areko::remesher::merge
