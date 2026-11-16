#include "pipeline.h"

#include "api/status.h"
#include "curves/curve_distribution.h"
#include "extract/extraction.h"
#include "field/border_lock.h"
#include "field/cross_field_solver.h"
#include "field/curvature_field.h"
#include "field/cut_graph_builder.h"
#include "field/part_edge_length.h"
#include "field/size_field_builder.h"
#include "input/input_builder.h"
#include "merge/result_assembly.h"
#include "param/chain_builder.h"
#include "param/solution_builder.h"
#include "param/system_builder.h"
#include "parts/part_builder.h"
#include "parts/part_curve_cut.h"
#include "prep/cleanup.h"
#include "prep/feature_lines.h"
#include "prep/preparation.h"
#include "prep/symmetry.h"
#include "prep/symmetry_detector.h"
#include "prep_cache.h"

#include <algorithm>
#include <array>
#include <condition_variable>
#include <exception>
#include <iterator>
#include <mutex>
#include <optional>
#include <thread>

namespace areko::remesher::run {

namespace {

// Per face 0 or 1 from the face sides of the mirror pairs: those of the cut,
// or, when the symmetry is to be found, those found on the mesh as it is now.
std::vector<int> mirrorSides(const prep::Prepared& prepared)
{
    std::optional<prep::SymmetryCache> found;
    const prep::SymmetryCache* cache = prepared.state.mirror ? &*prepared.state.mirror : nullptr;
    if (cache == nullptr && prepared.state.symmetry.detect) {
        found = prep::SymmetryDetector::run(prepared.mesh, nullptr);
        cache = found ? &*found : nullptr;
    }
    std::vector<int> sides;
    if (cache != nullptr) {
        const mesh::BitArray& faceSide = cache->faceSide();
        sides.resize(faceSide.size());
        for (std::size_t f = 0; f < faceSide.size(); ++f) {
            sides[f] = faceSide.test(f) ? 1 : 0;
        }
    }
    return sides;
}

// Guide curves of the host, in input order.
std::vector<curves::GuideCurve> hostCurves(const input::MeshInput& input)
{
    std::vector<curves::GuideCurve> out;
    for (std::size_t i = 0; i < input.curves().size(); ++i) {
        const input::MeshInput::Curve& source = input.curves()[i];
        curves::GuideCurve curve;
        curve.index = static_cast<int>(i);
        for (std::size_t k = 0; k + 2 < source.points.size(); k += 3) {
            curve.points.push_back({source.points[k], source.points[k + 1], source.points[k + 2]});
        }
        curve.kind = source.kind;
        curve.active = source.active;
        curve.width = source.width;
        curve.strength = source.strength / 100.0f;
        out.push_back(std::move(curve));
    }
    return out;
}

// The symmetry names exactly one mirror axis, without detection and without a kept half.
bool singleMirrorAxis(const settings::SymmetrySwitches& s)
{
    const int axes = static_cast<int>(s.axis[0]) + static_cast<int>(s.axis[1]) + static_cast<int>(s.axis[2]);
    return s.on && !s.detect && axes == 1;
}

// Number of threads for the parts: the setting, 0 meaning all cores.
int partThreads(int setting)
{
    if (setting > 0) {
        return setting;
    }
    const unsigned cores = std::thread::hardware_concurrency();
    return cores > 0 ? static_cast<int>(cores) : 1;
}

// The share done of one part, reported while the parts run one after the
// other; with several threads a part reports nothing (the share of finished
// parts is reported instead).
class PartReport {
public:
    explicit PartReport(Progress* progress) noexcept : progress_(progress) {}

    void operator()(float fraction) const
    {
        if (progress_ != nullptr) {
            progress_->report(fraction);
        }
    }

private:
    Progress* progress_;
};

// Runs `step(part, report)` for every part that is not skipped.
//  - One thread (or one part to compute): the parts run in order on this
//    thread, asking for a stop before each part.
//  - Several: up to `threads` threads take the parts in order, each part
//    computed by one thread alone (a step reads and writes only its own part
//    and reads the settings of the run, so every part comes out as in order).
//    A thread asks for a stop before it takes a part. After a failure no part
//    is taken any more, the parts already taken finish, and the failure of the
//    first part in order is passed on, as the run in order would have stopped
//    there. This thread only reports the share of finished parts, so the
//    progress callback stays on the thread of the run.
template <class Step>
void forEachComputedPart(std::vector<parts::Part>& list, Progress& progress, int threads, Step step)
{
    std::vector<parts::Part*> todo;
    for (parts::Part& part : list) {
        if (!part.skipped) {
            todo.push_back(&part);
        }
    }
    const std::size_t workers = std::min(static_cast<std::size_t>(threads > 1 ? threads : 1), todo.size());
    if (workers <= 1) {
        for (parts::Part& part : list) {
            progress.check();
            if (!part.skipped) {
                step(part, PartReport(&progress));
            }
        }
        return;
    }
    std::mutex mutex;
    std::condition_variable changed;
    std::size_t next = 0;
    std::size_t finished = 0;
    std::size_t stopped = 0;
    bool failed = false;
    std::vector<std::exception_ptr> errors(todo.size());
    const auto work = [&] {
        std::unique_lock lock(mutex);
        while (!failed && next < todo.size()) {
            const std::size_t i = next++;
            lock.unlock();
            std::exception_ptr error;
            try {
                progress.check();
                step(*todo[i], PartReport(nullptr));
            } catch (...) {
                error = std::current_exception();
            }
            lock.lock();
            if (error) {
                errors[i] = error;
                failed = true;
            } else {
                ++finished;
            }
            changed.notify_all();
        }
        ++stopped;
        changed.notify_all();
    };
    std::vector<std::thread> pool;
    pool.reserve(workers);
    std::exception_ptr startError;
    try {
        for (std::size_t w = 0; w < workers; ++w) {
            pool.emplace_back(work);
        }
    } catch (...) {
        startError = std::current_exception();
        const std::lock_guard lock(mutex);
        failed = true;
    }
    {
        std::unique_lock lock(mutex);
        std::size_t reported = 0;
        while (stopped < pool.size()) {
            changed.wait(lock, [&] { return stopped == pool.size() || finished != reported; });
            if (finished != reported) {
                reported = finished;
                const float share = static_cast<float>(reported) / static_cast<float>(todo.size());
                lock.unlock();
                progress.report(share);
                lock.lock();
            }
        }
    }
    for (std::thread& thread : pool) {
        thread.join();
    }
    for (const std::exception_ptr& error : errors) {
        if (error) {
            std::rethrow_exception(error);
        }
    }
    if (startError) {
        std::rethrow_exception(startError);
    }
}

// Stage data of the parts: one mesh (the triangle meshes of the parts, or
// with `rawMeshes` their raw meshes, a part without one adding nothing) and
// the part per face; for the stage "curves" also the cut chains as lines.
// Returns the first point of every part in the stage mesh.
std::vector<int> partStage(const std::vector<parts::Part>& list, bool withChains, bool rawMeshes, StageResult& result)
{
    const parts::PartBuilder::MeshOf triangles = [](const parts::Part& part) -> const mesh::Mesh* {
        return &part.mesh;
    };
    std::vector<int> pointStarts;
    result.mesh = parts::PartBuilder::merged(list, rawMeshes ? parts::PartBuilder::rawMeshOf : triangles,
                                             result.faceParts, pointStarts);
    result.lines.clear();
    if (!withChains) {
        return pointStarts;
    }
    for (std::size_t p = 0; p < list.size(); ++p) {
        const int offset = pointStarts[p];
        for (const curves::GuideCurve& curve : list[p].curves) {
            if (curve.cutVertices.size() < 2) {
                continue;
            }
            prep::NetLine line;
            for (const int v : curve.cutVertices) {
                line.vertices.push_back(v + offset);
            }
            line.points = curve.cutPoints;
            line.type = curve.kind;
            line.flags = curve.lineFlags;
            result.lines.push_back(std::move(line));
            result.lineKinds.push_back(curve.fromHost() ? ARM_LINE_GUIDE_CURVE : ARM_LINE_FEATURE);
        }
    }
    return pointStarts;
}

// One direction (xyz) per face of every part, in the face order of the merged
// stage mesh: `directionOf(part, t)` for a part where `has(part)` holds, zero
// for the faces of every other part.
template <class Has, class DirectionOf>
std::vector<float> faceDirections(const std::vector<parts::Part>& list, Has has, DirectionOf directionOf)
{
    std::vector<float> out;
    for (const parts::Part& part : list) {
        const std::size_t faces = static_cast<std::size_t>(part.mesh.faceCount());
        if (!has(part)) {
            out.insert(out.end(), 3 * faces, 0.0f);
            continue;
        }
        for (std::size_t t = 0; t < faces; ++t) {
            const mesh::Vec3f d = directionOf(part, t);
            out.insert(out.end(), {d.x, d.y, d.z});
        }
    }
    return out;
}

// Cross field of all parts in the face and point order of the merged stage
// mesh: one arm of the cross per face, and the singular points with the turn
// of the crosses around them in quarter turns.
void crossFieldStage(const std::vector<parts::Part>& list, const std::vector<int>& pointStarts, StageResult& result)
{
    result.crossDirections = faceDirections(
        list, [](const parts::Part& part) { return part.crossField.has_value(); },
        [](const parts::Part& part, std::size_t t) { return part.crossField->direction(t); });
    for (std::size_t p = 0; p < list.size(); ++p) {
        if (!list[p].crossField) {
            continue;
        }
        const field::Singularities& singular = list[p].crossField->singularities;
        for (std::size_t i = 0; i < singular.points.size(); ++i) {
            result.singularities.insert(result.singularities.end(),
                                        {singular.points[i] + pointStarts[p], singular.quarterTurns(i)});
        }
    }
}

// Cut edges of all parts as point pairs of the merged stage mesh.
std::vector<int> cutEdgeStage(const std::vector<parts::Part>& list, const std::vector<int>& pointStarts)
{
    std::vector<int> out;
    for (std::size_t p = 0; p < list.size(); ++p) {
        const parts::Part& part = list[p];
        if (!part.cutGraph) {
            continue;
        }
        const mesh::Adjacency& adjacency = part.mesh.adjacency();
        for (const int e : part.cutGraph->cutEdges) {
            const auto ends = adjacency.edgeVertices(e);
            out.insert(out.end(), {ends[0] + pointStarts[p], ends[1] + pointStarts[p]});
        }
    }
    return out;
}

// Size field of all parts: per face of the merged stage mesh its two target
// sizes (zero for a part without size field), and per part its target edge
// length, the smallest and largest size and the estimated number of quads.
void sizeFieldStage(const std::vector<parts::Part>& list, StageResult& result)
{
    for (const parts::Part& part : list) {
        const std::size_t faces = static_cast<std::size_t>(part.mesh.faceCount());
        if (!part.sizeField) {
            result.faceSizes.insert(result.faceSizes.end(), 2 * faces, 0.0f);
            result.partSizes.insert(result.partSizes.end(), {part.edgeLength, 0.0f, 0.0f, 0.0f});
            continue;
        }
        const field::SizeField& sizes = *part.sizeField;
        for (std::size_t t = 0; t < faces; ++t) {
            result.faceSizes.insert(result.faceSizes.end(), {sizes.sizeU[t], sizes.sizeV[t]});
        }
        result.partSizes.insert(result.partSizes.end(), {sizes.edgeLength, sizes.sizeMin, sizes.sizeMax,
                                                         static_cast<float>(sizes.estimatedQuads)});
    }
}

// Grid coordinates of all parts: per corner of the merged stage mesh the
// coordinates of its point of the cut copy (zero for a part without
// solution), and per part whether it has a solution.
void gridStage(const std::vector<parts::Part>& list, StageResult& result)
{
    for (const parts::Part& part : list) {
        const int corners = part.mesh.cornerCount();
        const bool solved = part.solution.has_value() && part.cutGraph.has_value();
        result.gridParts.push_back(solved ? 1 : 0);
        if (!solved) {
            result.gridCoordinates.insert(result.gridCoordinates.end(), 2 * static_cast<std::size_t>(corners), 0.0f);
            continue;
        }
        const mesh::Topology& copy = part.cutGraph->copy.topology();
        for (int c = 0; c < corners; ++c) {
            const std::array<float, 2>& uv = part.solution->uv[static_cast<std::size_t>(copy.vertex(c))];
            result.gridCoordinates.insert(result.gridCoordinates.end(), {uv[0], uv[1]});
        }
    }
}

// Curvature directions of all parts in the face order of the merged stage mesh.
std::vector<float> curvatureDirections(const std::vector<parts::Part>& list)
{
    return faceDirections(
        list, [](const parts::Part& part) { return part.curvature.has_value(); },
        [](const parts::Part& part, std::size_t t) { return part.curvature->direction[t]; });
}

// Stage data of the stages from "parts" to "quads": the merged stage mesh
// and what the reached stage shows on it.
void partStageData(const std::vector<parts::Part>& list, arm_stage until, StageResult& result)
{
    // The lines of the curves lie on the triangles of the parts; the raw mesh
    // has other points.
    const bool chains = until >= ARM_STAGE_CURVES && until < ARM_STAGE_RAW_MESH;
    const std::vector<int> pointStarts = partStage(list, chains, result.stage >= ARM_STAGE_RAW_MESH, result);
    if (result.stage == ARM_STAGE_CURVATURE) {
        result.faceDirections = curvatureDirections(list);
    } else if (result.stage == ARM_STAGE_CROSS_FIELD) {
        crossFieldStage(list, pointStarts, result);
    } else if (result.stage == ARM_STAGE_CUTS) {
        result.cutEdges = cutEdgeStage(list, pointStarts);
    } else if (result.stage == ARM_STAGE_SIZE_FIELD) {
        sizeFieldStage(list, result);
    } else if (result.stage == ARM_STAGE_PARAMETERIZATION) {
        result.cutEdges = cutEdgeStage(list, pointStarts);
        gridStage(list, result);
    }
}

} // namespace

StageResult Pipeline::run(arm_stage until)
{
    if (until < ARM_STAGE_INPUT || until > lastAvailable) {
        throw api::StatusError(ARM_ERROR_STAGE_NOT_AVAILABLE, "stage not available in this version");
    }
    if (!input_.hasMesh()) {
        throw api::StatusError(ARM_ERROR_NO_MESH, "no input mesh");
    }
    // One weight per stage that can be computed.
    static constexpr float weights[] = {1.0f, 9.0f, 9.0f, 1.0f, 2.0f, 1.0f, 2.0f, 4.0f, 6.0f, 2.0f, 3.0f, 12.0f,
                                        1.0f, 1.0f, 1.0f};
    static_assert(std::size(weights) == static_cast<std::size_t>(lastAvailable) + 1);
    progress_.plan(std::span(weights, static_cast<std::size_t>(until) + 1));
    progress_.beginPhase(0, ARM_STAGE_INPUT);

    StageResult result;
    result.stage = ARM_STAGE_INPUT;
    // The remembered preparation stands for the stages "input" and "prepared"
    // when its key matches; the input mesh has the faces of the input, so the
    // derived values are the same before it is built.
    std::optional<PrepKey> key;
    std::optional<prep::Prepared> remembered;
    if (cache_ != nullptr && until >= ARM_STAGE_PREPARED) {
        result.derived = settings::SettingsConversion::derive(settings_, input_.faceCount());
        key = PrepKey::of(input_, settings_, result.derived);
        remembered = cache_->find(*key);
        result.preparationReused = remembered.has_value();
    }
    if (!remembered) {
        result.mesh = input::InputBuilder::build(input_, settings_);
        progress_.check();
        result.derived = settings::SettingsConversion::derive(settings_, result.mesh.faceCount());
    }
    if (until >= ARM_STAGE_PREPARED) {
        progress_.beginPhase(1, ARM_STAGE_PREPARED);
        prep::Prepared prepared = remembered
                                      ? std::move(*remembered)
                                      : prep::Preparation::run(std::move(result.mesh), settings_, result.derived, progress_);
        remembered.reset();
        if (key && !result.preparationReused) {
            cache_->store(*key, prepared);
        }
        progress_.check();
        for (const input::MeshInput::Curve& curve : input_.curves()) {
            prepared.state.hasCuttingCurves = prepared.state.hasCuttingCurves || curves::cutsIntoMesh(curve.kind);
        }
        result.stage = ARM_STAGE_PREPARED;
        if (until >= ARM_STAGE_CLEANED) {
            progress_.check();
            progress_.beginPhase(2, ARM_STAGE_CLEANED);
            prep::Cleanup::run(prepared, progress_);
            result.stage = ARM_STAGE_CLEANED;
            if (until >= ARM_STAGE_MIRROR_HALF) {
                progress_.check();
                progress_.beginPhase(3, ARM_STAGE_MIRROR_HALF);
                prep::Symmetry::run(prepared, progress_);
                result.stage = ARM_STAGE_MIRROR_HALF;
                if (until == ARM_STAGE_MIRROR_HALF) {
                    result.mirrorSides = mirrorSides(prepared);
                }
                if (until >= ARM_STAGE_FEATURE_LINES) {
                    progress_.check();
                    progress_.beginPhase(4, ARM_STAGE_FEATURE_LINES);
                    prep::FeatureLines::run(prepared, progress_);
                    result.stage = ARM_STAGE_FEATURE_LINES;
                }
            }
            result.lines = prepared.state.lines.release();
        }
        if (until >= ARM_STAGE_PARTS) {
            progress_.check();
            progress_.beginPhase(5, ARM_STAGE_PARTS);
            const prep::PrepState& state = prepared.state;
            parts::PartBuilder::Split split = parts::PartBuilder::split(
                prepared.mesh, state.mirror ? &*state.mirror : nullptr, state.symmetry.on);
            if (split.sharedPoints) {
                throw api::StatusError(ARM_ERROR_NO_RESULT, "a point is shared by separate pieces of the mesh");
            }
            curves::CurveDistribution::run(prepared.mesh, hostCurves(input_), result.lines, split);
            parts::PartBuilder::assignBudgets(split.parts, prepared.mesh, state.quadBudget, state.targetEdgeLength);
            for (parts::Part& part : split.parts) {
                part.skipped = parts::PartBuilder::negligible(part, prepared.mesh);
            }
            result.stage = ARM_STAGE_PARTS;
            const int threads = partThreads(settings_.thread_count);
            if (until >= ARM_STAGE_CURVES) {
                progress_.beginPhase(6, ARM_STAGE_CURVES);
                const bool mirrorAxis = singleMirrorAxis(state.symmetry);
                forEachComputedPart(split.parts, progress_, threads, [&](parts::Part& part, const PartReport&) {
                    parts::PartCurveCut::run(part, mirrorAxis);
                });
                result.stage = ARM_STAGE_CURVES;
            }
            if (until >= ARM_STAGE_CURVATURE) {
                progress_.beginPhase(7, ARM_STAGE_CURVATURE);
                field::CurvatureField::Options options;
                options.symmetry = state.symmetry.on;
                options.halfComputed = state.symmetry.curvatureSymmetry;
                options.followBorders = settings_.follow_borders != 0;
                options.severalParts = split.parts.size() >= 2;
                field::BorderLock::Options lock = field::BorderLock::optionsFor(state.symmetry);
                if (lock.halfMode) {
                    lock.planeMax = prepared.mesh.bounds().max;
                }
                forEachComputedPart(split.parts, progress_, threads, [&](parts::Part& part, const PartReport&) {
                    field::BorderLock::assign(part, lock);
                    field::PartEdgeLength::assign(part, state.quadCount, split.parts.size(),
                                                  result.derived.borderDensityInverse);
                    field::CurvatureField::compute(part, options, progress_);
                });
                result.stage = ARM_STAGE_CURVATURE;
            }
            if (until >= ARM_STAGE_CROSS_FIELD) {
                progress_.beginPhase(8, ARM_STAGE_CROSS_FIELD);
                field::CrossFieldSolver::Options options;
                options.symmetry = state.symmetry.on;
                forEachComputedPart(split.parts, progress_, threads, [&](parts::Part& part, const PartReport&) {
                    field::CrossFieldSolver::compute(part, options, progress_);
                });
                result.stage = ARM_STAGE_CROSS_FIELD;
            }
            if (until >= ARM_STAGE_CUTS) {
                progress_.beginPhase(9, ARM_STAGE_CUTS);
                field::CutGraphBuilder::Options options;
                options.symmetry = state.symmetry.on;
                options.followBorders = settings_.follow_borders != 0;
                forEachComputedPart(split.parts, progress_, threads, [&](parts::Part& part, const PartReport&) {
                    // A part without cross field ends here: it gets no cut graph.
                    if (part.crossField) {
                        field::CutGraphBuilder::compute(part, options, progress_);
                    }
                });
                result.stage = ARM_STAGE_CUTS;
            }
            if (until >= ARM_STAGE_SIZE_FIELD) {
                progress_.beginPhase(10, ARM_STAGE_SIZE_FIELD);
                field::SizeFieldBuilder::Options options;
                options.exactCount = settings_.exact_count != 0;
                options.pointsPerCircle = result.derived.pointsPerCircle;
                options.densityRatio = result.derived.densityRatio;
                options.maxQuadRatio = settings_.max_quad_ratio;
                options.borderDensityInverse = result.derived.borderDensityInverse;
                options.colorDensity = settings_.color_density != 0;
                options.valueDensity = settings_.attribute_density != 0;
                forEachComputedPart(split.parts, progress_, threads, [&](parts::Part& part, const PartReport&) {
                    field::PartEdgeLength::assign(part, state.quadCount, split.parts.size(),
                                                  options.borderDensityInverse);
                    field::SizeFieldBuilder::compute(part, options, progress_);
                });
                result.stage = ARM_STAGE_SIZE_FIELD;
            }
            if (until >= ARM_STAGE_PARAMETERIZATION) {
                progress_.beginPhase(11, ARM_STAGE_PARAMETERIZATION);
                param::ChainBuilder::Options options;
                options.symmetry = state.symmetry.on;
                // Share of a part's time that is over after its chains and
                // after its system; solving takes the rest.
                static constexpr float afterChains = 0.02f;
                static constexpr float afterSystem = 0.05f;
                const float partCount = static_cast<float>(split.parts.size());
                forEachComputedPart(split.parts, progress_, threads, [&](parts::Part& part, const PartReport& report) {
                    // A part without cut graph has no field to follow: it gets no chains.
                    if (part.cutGraph) {
                        const float before = static_cast<float>(part.index);
                        param::ChainBuilder::compute(part, options, progress_);
                        report((before + afterChains) / partCount);
                        param::SystemBuilder::Options system;
                        system.symmetry = state.symmetry.on;
                        param::SystemBuilder::compute(part, system, progress_);
                        report((before + afterSystem) / partCount);
                        param::SolutionBuilder::Options solution;
                        solution.symmetry = state.symmetry.on;
                        solution.freezeBorders = settings_.freeze_borders != 0;
                        solution.report = [&, before](float solved) {
                            report((before + afterSystem + (1.0f - afterSystem) * solved) / partCount);
                        };
                        param::SolutionBuilder::compute(part, solution, progress_);
                    }
                });
                result.stage = ARM_STAGE_PARAMETERIZATION;
            }
            extract::Extraction::Options extraction;
            extraction.followBorders = settings_.follow_borders != 0;
            extraction.keepColors = settings_.keep_colors != 0;
            extraction.symmetry = state.symmetry.on;
            const float partCount = static_cast<float>(split.parts.size());
            if (until >= ARM_STAGE_RAW_MESH) {
                progress_.beginPhase(12, ARM_STAGE_RAW_MESH);
                forEachComputedPart(split.parts, progress_, threads, [&](parts::Part& part, const PartReport& report) {
                    extract::Extraction::rawMesh(part, extraction, progress_);
                    report(static_cast<float>(part.index + 1) / partCount);
                });
                result.stage = ARM_STAGE_RAW_MESH;
            }
            if (until >= ARM_STAGE_QUADS) {
                progress_.beginPhase(13, ARM_STAGE_QUADS);
                forEachComputedPart(split.parts, progress_, threads, [&](parts::Part& part, const PartReport& report) {
                    extract::Extraction::quads(part, extraction);
                    report(static_cast<float>(part.index + 1) / partCount);
                });
                result.stage = ARM_STAGE_QUADS;
            }
            // Stage "result": the features of the triangle mesh follow onto
            // the quad mesh of every part, the quad meshes are finished, then
            // they become the one result mesh.
            if (until >= ARM_STAGE_RESULT) {
                progress_.beginPhase(14, ARM_STAGE_RESULT);
                forEachComputedPart(split.parts, progress_, threads, [&](parts::Part& part, const PartReport& report) {
                    extract::Extraction::finish(part, extraction, progress_);
                    report(static_cast<float>(part.index + 1) / partCount);
                });
                progress_.check();
                merge::ResultAssembly::Options assembly;
                assembly.keepColors = settings_.keep_colors != 0;
                assembly.splitAtCreases = settings_.split_at_creases != 0;
                assembly.cleanOutput = settings_.clean_output != 0;
                result.mesh = merge::ResultAssembly::run(split.parts, prepared, assembly);
                result.stage = ARM_STAGE_RESULT;
            }
            for (const parts::Part& part : split.parts) {
                result.partValues.insert(result.partValues.end(), {part.quadBudget, part.skipped ? 1 : 0});
            }
            if (result.stage == ARM_STAGE_RESULT) {
                // The result has no lines and, mirrored and split, no part per face.
                result.lines.clear();
            } else {
                partStageData(split.parts, until, result);
            }
        } else {
            result.mesh = std::move(prepared.mesh);
        }
    }
    progress_.finish();
    return result;
}

} // namespace areko::remesher::run
