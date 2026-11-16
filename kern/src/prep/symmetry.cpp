#include "symmetry.h"

#include "mesh_repair.h"
#include "symmetry_cut.h"
#include "triangulator.h"

#include <utility>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

void Symmetry::run(Prepared& prepared, run::Progress& progress)
{
    Mesh& mesh = prepared.mesh;
    PrepState& state = prepared.state;
    bool mirrorTriangulation = false;
    if (state.symmetry.on) {
        progress.check();
        if (state.hasEdgeAttribute) {
            state.lines.toSharpEdges(mesh);
        }
        const bool* axes = state.symmetry.axis;
        const int k = static_cast<int>(axes[0]) + static_cast<int>(axes[1]) + static_cast<int>(axes[2]);
        const bool halfOnly = k > 1 || state.hasCuttingCurves;
        for (int a = 0; a < 3; ++a) {
            if (!axes[a]) {
                continue;
            }
            SymmetryCut::Job job;
            job.axis = a;
            job.splitEdgesFirst = true;
            job.cutOnly = halfOnly;
            job.local = state.symmetryLocal;
            job.usePresetPlane = true;
            job.planePos = 0.0f;
            job.keepPointOrigin = true;
            job.keepFaceOrigin = PrepState::useGroupBorders;
            job.keepHardEdges = PrepState::useEdgeAttribute;
            SymmetryCut::Result result = SymmetryCut::run(mesh, job);
            state.symmetryPlane[static_cast<std::size_t>(a)] = result.planePos;
            if (result.cache) {
                state.mirror = std::move(result.cache);
            }
        }
        if (halfOnly) {
            state.symmetry.halfMode = true;
            state.quadBudget = state.quadCount / (1 << k);
        }
        mirrorTriangulation = k != 0;
        state.lines.extractFromSharpEdges(mesh);
    }
    progress.check();

    if (mesh.topology().mode() == FaceMode::Polygons) {
        const MirrorLists lists = state.mirror ? state.mirror->mirrorLists() : MirrorLists{};
        Triangulator::Options options;
        options.mirror = mirrorTriangulation ? &lists : nullptr;
        mesh = Triangulator::run(std::move(mesh), options);
        if (!mesh.adjacency().nonManifoldEdges().empty()) {
            (void)MeshRepair::run(mesh);
        }
        // The mirror pairs pass over to the triangles as they are.
    }
}

} // namespace areko::remesher::prep
