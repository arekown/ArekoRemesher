#include "cleanup.h"

#include "decimator.h"
#include "degenerate_groups.h"
#include "denoiser.h"
#include "freeze_marks.h"
#include "line_refresh.h"
#include "mesh_cleaner.h"
#include "refiner.h"
#include "small_region_merge.h"

#include "mesh/sharp_edges.h"

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

void Cleanup::run(Prepared& prepared, run::Progress& progress)
{
    Mesh& mesh = prepared.mesh;
    PrepState& state = prepared.state;
    constexpr int workResolution = PrepState::workResolution;

    if (state.denoiseStrength > 0.0f && mesh.pointCount() > workResolution) {
        Denoiser::run(mesh, state);
    }
    mesh.maps().remove(NetLines::originMap, ElementKind::Point);
    state.lines.extractFromSharpEdges(mesh);
    MeshCleaner::run(mesh, 0, state);
    state.lines.restore(mesh, NetLines::Restore::Direct, nullptr);
    state.lines.toSharpEdges(mesh);
    DegenerateGroups::run(mesh, state.referenceLength);
    SmallRegionMerge::run(mesh, state);
    LineRefresh::run(mesh, state);
    progress.check();

    if (mesh.pointCount() > workResolution) {
        state.lines.toSharpEdges(mesh);
        mesh = Decimator(state.decimationLength, workResolution, state.hasEdgeAttribute, progress).run(std::move(mesh));
        if (state.hasEdgeAttribute) {
            // Both forms of the sharp edges describe the same edges again.
            SharpEdges::toEdgeForm(mesh);
        }
    }
    state.surfaceArea = -1.0f;
    state.referenceLength = 0.0f;
    if (PrepState::useGroupBorders || state.hasEdgeAttribute) {
        state.lines.extractFromSharpEdges(mesh);
    }
    if (state.freezeBorders) {
        FreezeMarks::mark(mesh, true);
    }
    progress.check();

    MeshCleaner::run(mesh, 1, state);
    progress.check();

    Refiner::run(mesh, state, state.surfaceArea, state.refinePasses);
    progress.check();
}

} // namespace areko::remesher::prep
