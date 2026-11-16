#pragma once

#include "net_lines.h"
#include "symmetry_cache.h"

#include "settings/settings.h"

#include <array>
#include <optional>

namespace areko::remesher::prep {

// State of the preparation that its steps read and write, one per run.
struct PrepState {
    static constexpr int workResolution = 50000;
    // Face group borders and the sharp edge attribute are always used.
    static constexpr bool useGroupBorders = true;
    static constexpr bool useEdgeAttribute = true;

    int quadCount = 0;                 // target number of quads
    float targetEdgeLength = -1.0f;
    bool freezeBorders = false;
    float creaseAngle = 0.5235988f;   // radians
    bool autoDetectCreases = true;     // detected crease lines join the lines
    float denoiseStrength = 0.5f;      // only > 0 / = 0 matters
    bool hasEdgeAttribute = false;     // the input brought sharp edges
    float surfaceArea = -1.0f;         // area cache of the reference length
    float referenceLength = 0.0f;      // L
    float decimationLength = 0.0f;     // Ld
    int repairCount = -1;
    int refineVertexBudget = -1;       // < 1: 1.6 times the work resolution
    int refinePasses = 4;
    int refineInsertedVertexCount = -1;
    NetLines lines;

    // Symmetry: switches (the half mode is also set by the cut of several
    // axes), plane through the bounding box middle, plane position per axis,
    // quad budget, guiding curves that cut into the mesh, mirror pairs.
    settings::SymmetrySwitches symmetry;
    bool symmetryLocal = true;
    std::array<float, 3> symmetryPlane{0.0f, 0.0f, 0.0f};
    int quadBudget = 0;
    bool hasCuttingCurves = false;
    std::optional<SymmetryCache> mirror;

    // The stored area, or the area of the mesh while none is stored (below 0).
    [[nodiscard]] float knownArea(const mesh::Mesh& mesh) const
    {
        return surfaceArea < 0.0f ? mesh.triangleGeometry().totalArea() : surfaceArea;
    }
};

} // namespace areko::remesher::prep
