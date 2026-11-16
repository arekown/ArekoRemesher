#pragma once

#include "areko_remesher.h"
#include "input/mesh_input.h"
#include "mesh/mesh.h"
#include "prep/net_lines.h"
#include "progress.h"
#include "settings/settings.h"

namespace areko::remesher::run {

struct StageResult {
    arm_stage stage = ARM_STAGE_INPUT;
    // The stages "input" and "prepared" were skipped: their result came from
    // the remembered preparation.
    bool preparationReused = false;
    mesh::Mesh mesh;
    settings::DerivedSettings derived;
    // Lines of the preparation (from the stage "cleaned" on).
    std::vector<prep::NetLine> lines;
    // Per face 0 (computed half) or 1 (mirrored half) when the stage "mirror
    // half" is shown and mirrored the mesh or found its symmetry, else empty.
    std::vector<int> mirrorSides;
    // From the stage "parts" on: part number per face, and per part its quad
    // budget and 1 when it is skipped (pairs); kind per line (arm_line_kind).
    std::vector<int> faceParts;
    std::vector<int> partValues;
    std::vector<int> lineKinds;
    // Stage "curvature": per face the direction of least curvature (xyz),
    // zero where none was found.
    std::vector<float> faceDirections;
    // Stage "cross field": per face one arm of its cross (xyz, zero where the
    // part has no cross field), and per singular point its number in the stage
    // mesh and the turn of the crosses around it in quarter turns (pairs).
    std::vector<float> crossDirections;
    std::vector<int> singularities;
    // Stage "cuts": the edges along which the parts are cut open, as point
    // pairs of the stage mesh.
    std::vector<int> cutEdges;
    // Stage "size field": per face its target size along the two directions
    // of its cross (pairs), and per part its target edge length, the smallest
    // and largest size and the estimated number of quads (four values).
    std::vector<float> faceSizes;
    std::vector<float> partSizes;
    // Stage "parameterization": per corner its two grid coordinates (pairs;
    // zero for a part without solution), and per part 1 when it has a solution.
    std::vector<float> gridCoordinates;
    std::vector<int> gridParts;
    // Stages "raw mesh" and "quads": the mesh is the polygon mesh of all
    // parts (unused points included), with the part number per face.
    // Stage "result": the finished mesh with its maps (materials, face
    // groups, smoothing groups, sharp edges, colours), without part per face.
};

class PrepCache;

// Runs the stages in order up to the requested one. With a cache, a run that
// reaches the stage "prepared" takes its result from the cache when the key
// matches, and otherwise leaves its own result there.
class Pipeline {
public:
    static constexpr arm_stage lastAvailable = ARM_STAGE_RESULT;

    Pipeline(const input::MeshInput& input, const arm_settings& settings, Progress& progress,
             PrepCache* cache = nullptr) noexcept
        : input_(input), settings_(settings), progress_(progress), cache_(cache)
    {
    }

    [[nodiscard]] StageResult run(arm_stage until);

private:
    const input::MeshInput& input_;
    const arm_settings& settings_;
    Progress& progress_;
    PrepCache* cache_;
};

} // namespace areko::remesher::run
