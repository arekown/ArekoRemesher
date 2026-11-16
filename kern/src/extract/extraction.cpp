#include "extraction.h"

#include "border_classifier.h"
#include "doublets.h"
#include "face_closing.h"
#include "feature_transfer.h"
#include "hole_closing.h"
#include "isoline_tracer.h"
#include "ngon_optimizer.h"
#include "quad_mesh_builder.h"
#include "final_cleanup.h"
#include "final_relaxation.h"
#include "frozen_border_reset.h"
#include "vertex_color_transfer.h"

#include "mesh/degenerate_cleanup.h"
#include "parts/part.h"

#include <optional>

namespace areko::remesher::extract {

void Extraction::rawMesh(parts::Part& part, const Options& options, const run::Progress& progress)
{
    part.raw.reset();
    const std::vector<mesh::FreezeMark>* marks = part.frozenBorderMarks();
    FaceClosing closer(marks);
    IsolineTracer::trace(part, closer, progress);
    if (!part.trace) {
        return;
    }
    progress.check();
    RawMesh raw = QuadMeshBuilder::build(part, *part.trace);
    progress.check();
    BorderClassifier::run(raw, part);
    mesh::DegenerateCleanup::run(raw.mesh, {});
    Doublets::run(raw.mesh);
    progress.check();
    std::optional<FrozenBorderInput> frozen;
    if (marks != nullptr) {
        frozen.emplace(FrozenBorderInput{*marks, raw.originalToOutput, &part.lockedBorderEdges});
    }
    const HoleClosing::Options holes{.rebuildPolygons = true, .followBorders = options.followBorders};
    HoleClosing(raw.mesh, raw.points, part.mesh, part.quadBudget, frozen ? &*frozen : nullptr, holes).run();
    part.raw = std::move(raw);
}

void Extraction::quads(parts::Part& part, const Options& options)
{
    if (!part.raw) {
        return;
    }
    RawMesh& raw = *part.raw;
    NgonOptimizer::run(raw.mesh, NgonSettings::forPart(part.frozenBorderMarks() != nullptr, options.followBorders));
    // The point data follow the new points.
    raw.points.resize(static_cast<std::size_t>(raw.mesh.pointCount()), PointSource{0, 0}, 0.0f);
}

void Extraction::finish(parts::Part& part, const Options& options, const run::Progress& progress)
{
    if (!part.raw) {
        return;
    }
    features(part, options);
    RawMesh& raw = *part.raw;
    FinalRelaxation::run(part, options.symmetry);
    progress.check();
    if (const std::vector<mesh::FreezeMark>* marks = part.frozenBorderMarks()) {
        FrozenBorderReset::run(raw.mesh, part.mesh, *marks, raw.originalToOutput);
    }
    progress.check();
    if (options.keepColors) {
        VertexColorTransfer::run(raw.mesh, raw.points.source, part.mesh);
    }
    FinalCleanup::run(raw);
}

void Extraction::features(parts::Part& part, const Options& options)
{
    FeatureTransfer::run(part, options.cutAlongCurves);
}

} // namespace areko::remesher::extract
