#include "feature_lines.h"

#include "denoiser.h"
#include "feature_chains.h"

#include "mesh/attribute_map.h"
#include "mesh/attribute_maps.h"

#include <memory>
#include <optional>
#include <string>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

void FeatureLines::run(Prepared& prepared, run::Progress& progress)
{
    Mesh& mesh = prepared.mesh;
    PrepState& state = prepared.state;

    const std::optional<NetLines::FinalLines> finalLines = state.lines.buildFinal(
        mesh, {.sharpEdges = PrepState::useEdgeAttribute, .groupBorders = PrepState::useGroupBorders});
    if (finalLines && finalLines->linePoints.size() != 0) {
        mesh.maps().set(std::make_unique<BoolMap>(std::string(NetLines::pointMap), ElementKind::Point,
                                                  finalLines->linePoints));
    }
    progress.check();

    if (state.denoiseStrength > 0.0f || state.autoDetectCreases) {
        state.surfaceArea = state.knownArea(mesh);
        Denoiser::Options options;
        options.mode = Denoiser::Mode::AfterDetection;
        options.chains.creaseAngle = state.creaseAngle;
        options.chains.lengths = FeatureChains::referenceLengths(state.surfaceArea, state.quadCount);
        if (PrepState::useGroupBorders || state.hasEdgeAttribute) {
            if (const auto* points = mesh.maps().findAs<BoolMap>(NetLines::pointMap, ElementKind::Point)) {
                options.fixedPoints = &points->bits();
            }
            // Without final lines there is no protected edge.
            options.chains.protectedEdges = finalLines ? &finalLines->protectedEdges : nullptr;
        }
        const FeatureChains::Result chains = Denoiser::run(mesh, options);
        progress.check();
        if (state.autoDetectCreases) {
            state.lines.markDetected(FeatureChains::CreaseLineCheck(mesh));
            state.lines.appendChains(mesh, chains.points, chains.classes, 4);
        }
    }
    if (!state.autoDetectCreases) {
        state.lines.detectedToProtection();
    }
    progress.check();
}

} // namespace areko::remesher::prep
