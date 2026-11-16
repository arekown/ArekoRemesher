#include "preparation.h"

#include "decimator.h"
#include "freeze_marks.h"
#include "import_attributes.h"
#include "reference_length.h"
#include "triangulator.h"

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

Prepared Preparation::run(Mesh&& input, const arm_settings& settings, const settings::DerivedSettings& derived,
                          run::Progress& progress)
{
    Prepared out;
    PrepState& state = out.state;
    out.groups = ImportAttributes::apply(input, settings);
    state.quadCount = derived.quadCount;
    state.targetEdgeLength = derived.edgeLength;
    state.freezeBorders = settings.freeze_borders != 0;
    state.creaseAngle = derived.creaseAngleRadians;
    state.autoDetectCreases = settings.detect_creases != 0;
    state.denoiseStrength = derived.denoiseStrength;
    state.hasEdgeAttribute = input.maps().find(map_names::sharpEdges) != nullptr;
    state.lines.reset(state.hasEdgeAttribute);
    state.symmetry = derived.symmetry;
    state.symmetryLocal = settings.mirror_at_center != 0;
    state.quadBudget = derived.quadCount;
    const bool freeze = state.freezeBorders;
    constexpr int workResolution = PrepState::workResolution;
    progress.check();

    Mesh triangles = Triangulator::run(std::move(input), Triangulator::Options{});
    if (freeze) {
        FreezeMarks::mark(triangles, true);
    }
    state.referenceLength =
        ReferenceLength::compute(triangles, derived.edgeLength, derived.quadCount, freeze, &state.surfaceArea);
    state.decimationLength = ReferenceLength::decimationLength(triangles, state.referenceLength, freeze, workResolution);
    progress.check();

    const int threshold = truncToInt(static_cast<float>(workResolution) * 1.7f);
    if (triangles.pointCount() > threshold) {
        const int target = truncToInt(static_cast<double>(workResolution) * 1.5);
        out.mesh = Decimator(state.decimationLength, target, state.hasEdgeAttribute, progress).run(std::move(triangles));
    } else {
        out.mesh = std::move(triangles);
    }
    return out;
}

} // namespace areko::remesher::prep
