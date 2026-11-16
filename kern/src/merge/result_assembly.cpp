#include "result_assembly.h"

#include "crease_split.h"
#include "group_restore.h"

#include "api/status.h"
#include "parts/part_builder.h"
#include "prep/edge_dissolver.h"
#include "prep/mesh_repair.h"
#include "prep/symmetry_cut.h"

#include <cstddef>
#include <set>
#include <string_view>
#include <vector>

namespace areko::remesher::merge {

using namespace areko::remesher::mesh;

namespace {

// New face groups for the parts whose face group map is empty: each takes
// the smallest value from 0 that no part has and that was not given before.
class GroupIds {
public:
    explicit GroupIds(std::vector<parts::Part>& parts)
    {
        for (parts::Part& part : parts) {
            if (!part.raw) {
                continue;
            }
            if (const TypedMap<int>* groups = part.raw->mesh.maps().faceGroups()) {
                used_.insert(groups->values().begin(), groups->values().end());
            }
        }
    }

    void assign(Mesh& mesh)
    {
        TypedMap<int>* groups = mesh.maps().faceGroups();
        if (groups == nullptr || !groups->values().empty()) {
            return;
        }
        int id = 0;
        while (used_.contains(id)) {
            ++id;
        }
        used_.insert(id);
        groups->values().assign(static_cast<std::size_t>(mesh.faceCount()), id);
    }

private:
    std::set<int> used_;
};

// Cuts the result at the plane of every mirror axis and mirrors the kept half.
void mirrorAtEnd(Mesh& mesh, const prep::PrepState& state, bool keepColors)
{
    for (int a = 0; a < 3; ++a) {
        if (!state.symmetry.axis[a]) {
            continue;
        }
        prep::SymmetryCut::Job job;
        job.axis = a;
        job.local = false;
        job.usePresetPlane = true;
        job.planePos = state.symmetryPlane[static_cast<std::size_t>(a)];
        job.buildCache = false;
        job.keepPointOrigin = keepColors;
        job.keepFaceOrigin = true;
        job.keepHardEdges = true;
        (void)prep::SymmetryCut::run(mesh, job);
    }
}

// Points between exactly two edges are dissolved, then non-manifold points split.
void cleanOutput(Mesh& mesh)
{
    (void)prep::EdgeDissolver::run(mesh);
    (void)prep::MeshRepair::splitNonManifoldPoints(mesh, nullptr);
}

// Materials, face groups and smoothing groups reach the host only with one
// value per face; a map of another length is dropped.
void dropFaceValuesOfOtherLength(Mesh& mesh)
{
    const auto faces = static_cast<std::size_t>(mesh.faceCount());
    for (const std::string_view name : {map_names::materials, map_names::faceGroups, map_names::smoothingGroups}) {
        const AttributeMap* map = mesh.maps().find(name, ElementKind::Face);
        if (map != nullptr && map->size() != faces) {
            mesh.maps().remove(name, ElementKind::Face);
        }
    }
}

} // namespace

Mesh ResultAssembly::run(std::vector<parts::Part>& parts, const prep::Prepared& prepared, const Options& options)
{
    if (parts.size() >= 2) {
        GroupIds ids(parts);
        for (parts::Part& part : parts) {
            if (part.raw) {
                ids.assign(part.raw->mesh);
            }
        }
    }
    std::vector<int> faceParts;
    std::vector<int> pointStarts;
    Mesh result = parts::PartBuilder::merged(parts, parts::PartBuilder::rawMeshOf, faceParts, pointStarts,
                                             {.result = true, .unite = true});
    result.toPolygonMode();
    const prep::PrepState& state = prepared.state;
    if (state.symmetry.on && state.symmetry.mirrorAtEnd) {
        mirrorAtEnd(result, state, options.keepColors);
    }
    if (prepared.groups.restore) {
        GroupRestore::run(result, prepared.groups);
    }
    if (options.splitAtCreases) {
        (void)CreaseSplit::run(result);
    }
    if (options.cleanOutput) {
        cleanOutput(result);
    }
    // The point colours stay with the result.
    if (result.faceCount() == 0) {
        throw api::StatusError(ARM_ERROR_NO_RESULT, "the remeshing left no faces");
    }
    dropFaceValuesOfOtherLength(result);
    return result;
}

} // namespace areko::remesher::merge
