#pragma once

#include "mesh/mesh.h"

#include <cstdint>
#include <vector>

namespace areko::remesher::prep {

// Values of the host per face group, so that the result can carry the host's
// materials, face groups and smoothing groups again. Index = face group value.
// A vector is empty when the host gave no such map; a fresh table belongs to
// every run.
struct GroupRestoreTable {
    std::vector<int> groupFaceGroup;
    std::vector<int> groupMaterial;
    std::vector<std::uint32_t> groupSmoothing;
    bool materialsAreGroups = false; // face groups are a copy of the materials
    bool restore = false;            // the result maps face groups back
};

// Turns materials, face groups and smoothing groups into face groups (regions)
// and sharp edges:
//  - neither materials nor smoothing groups: nothing changes,
//  - materials only (no face groups, no smoothing groups): face groups become a
//    copy of the materials,
//  - otherwise regions are the face components separated by edges between
//    different face groups, different materials and smoothing groups without a
//    common bit; a smoothing border inside one region becomes a sharp edge
//    (only edges with two corners that are no face group / material border).
//    Face groups become the region numbers.
class GroupUnion {
public:
    [[nodiscard]] static GroupRestoreTable apply(mesh::Mesh& mesh);
};

} // namespace areko::remesher::prep
