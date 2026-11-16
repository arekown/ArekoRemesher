#pragma once

#include "areko_remesher.h"
#include "group_union.h"
#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Turns the host attributes of the input into face groups and sharp edges, in
// this order: maps of switched-off sources are removed (materials, face
// groups, smoothing groups, sharp edges), group union (always), merging of CAD
// faces, normal splits, welding of open seams.
class ImportAttributes {
public:
    [[nodiscard]] static GroupRestoreTable apply(mesh::Mesh& mesh, const arm_settings& settings);
};

} // namespace areko::remesher::prep
