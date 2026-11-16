#include "mesh_cleaner.h"

#include "close_point_merge.h"
#include "duplicate_faces.h"
#include "edge_dissolver.h"
#include "hole_filler.h"
#include "mesh_repair.h"
#include "reference_length.h"
#include "triangulator.h"

#include <algorithm>
#include <bit>
#include <cstdint>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr int minimumQuads = 3000;
// Lower limit of the squared near distance: one step below 1e-18 in float.
const float minimumNearSquared = std::bit_cast<float>(std::uint32_t{0x219392ee});

} // namespace

void MeshCleaner::run(Mesh& mesh, int mode, PrepState& state)
{
    if (mode == 0) {
        (void)DuplicateFaces::remove(mesh);
    }

    const float boxDiagonal = diagonal(mesh.bounds());
    const float nearDistance = static_cast<float>(static_cast<double>(boxDiagonal) * 1e-8);
    const float nearSquared =
        static_cast<double>(nearDistance) < 1e-9 ? minimumNearSquared : nearDistance * nearDistance;
    const float smallDistance = static_cast<float>(static_cast<double>(boxDiagonal) * 5e-7);
    (void)ClosePointMerge::apply(mesh, nearSquared, smallDistance * smallDistance);

    const int firstRepair = MeshRepair::run(mesh);

    const float reference = ReferenceLength::compute(mesh, state.targetEdgeLength,
                                                     std::max(state.quadCount, minimumQuads), false, &state.surfaceArea);
    HoleFiller::Options holes;
    holes.maxPerimeter = 0.25f * reference;
    holes.maxArea = (reference * reference) / 50.0f;
    (void)HoleFiller::fill(mesh, holes);

    const int secondRepair = MeshRepair::run(mesh);

    const int dissolved = EdgeDissolver::run(mesh);
    if (dissolved > 0) {
        (void)MeshRepair::removeUnusedPoints(mesh);
    }
    state.repairCount = dissolved + firstRepair + secondRepair;

    if (mesh.topology().mode() != FaceMode::Triangles) {
        Triangulator::Options triangulation;
        triangulation.angleRule = false;
        mesh = Triangulator::run(std::move(mesh), triangulation);
    }
}

} // namespace areko::remesher::prep
