#pragma once

#include "prep_state.h"

#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Refinement of coarse working meshes (triangle meshes; the cleaning leaves
// one): flat triangle runs first, then up to `passes` edge splitting passes
// with the target lengths 2 Lr, Lr, Lr, 2 Lr (Lr from the target edge length
// or from `area` with the quad count clamped to 3000 .. 12000; with frozen
// borders bounded by the border length), a budget of 1.6 times the work
// resolution unless the state names one, the hard edge set renewed after every
// pass with new points (lines over two edges), a flat triangle run after the
// passes that need one, and a last run at the end. `area` is a copy: a newly
// computed area stays local. Writes the inserted point count of the state.
class Refiner {
public:
    static void run(mesh::Mesh& mesh, PrepState& state, float area, int passes);
};

} // namespace areko::remesher::prep
