#pragma once

#include "areko_remesher.h"
#include "group_union.h"
#include "prep_state.h"
#include "mesh/mesh.h"
#include "run/progress.h"
#include "settings/settings.h"

namespace areko::remesher::prep {

// Result of the stage "prepared": the triangle mesh after the first
// decimation, the state of the preparation and what later stages need from
// the import.
struct Prepared {
    mesh::Mesh mesh;
    GroupRestoreTable groups;
    PrepState state;
};

// First part of the preparation: import attributes, triangulation, freeze
// marks, reference lengths and, for meshes with more than 1.7 times the work
// resolution of points, the decimation to 1.5 times the work resolution.
class Preparation {
public:
    [[nodiscard]] static Prepared run(mesh::Mesh&& input, const arm_settings& settings,
                                      const settings::DerivedSettings& derived, run::Progress& progress);
};

} // namespace areko::remesher::prep
