#include "refiner.h"

#include "edge_splitter.h"
#include "flat_triangle_check.h"
#include "freeze_marks.h"
#include "reference_length.h"

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr int firstRuns = 7;
constexpr int defaultPasses = 4;
constexpr float budgetFactor = 1.6f;

// Lr: L0 from the target edge length or sqrt(area / clamped quad count); with
// frozen borders the border length Lb when 0.5 L0 < Lb <= L0, L0 when Lb > L0
// or Lb = 0, otherwise 0.5 L0.
float targetLength(const Mesh& mesh, const PrepState& state, float area)
{
    const float l0 = ReferenceLength::compute(mesh, state.targetEdgeLength,
                                              ReferenceLength::clampedQuadCount(state.quadCount), false, &area);
    if (!state.freezeBorders) {
        return l0;
    }
    const auto* marks = FreezeMarks::find(mesh);
    const float lb = ReferenceLength::borderLength(mesh, marks != nullptr ? &marks->values() : nullptr);
    if (lb > l0) {
        return l0;
    }
    if (lb > 0.5f * l0) {
        return lb;
    }
    if (lb == 0.0f) {
        return l0;
    }
    return 0.5f * l0;
}

// Runs the flat triangle check once; flipped edges renew the hard edge set.
FlatTriangleBlock flatRun(Mesh& mesh, PrepState& state, BitArray& hard, int mode, int round, int begin, int end)
{
    FlatTriangleBlock block;
    block.mode = mode;
    block.round = round;
    block.requiredBegin = begin;
    block.requiredEnd = end;
    FlatTriangleCheck::run(mesh, hard, block);
    if (block.count > 0 && mode != 2) {
        hard = state.lines.collectEdges(mesh, PrepState::useGroupBorders, state.hasEdgeAttribute, true);
    }
    return block;
}

} // namespace

void Refiner::run(Mesh& mesh, PrepState& state, float area, int passes)
{
    if (mesh.topology().mode() != FaceMode::Triangles) {
        return;
    }
    BitArray hard = state.lines.collectEdges(mesh, true, true, true);
    if (passes == 0) {
        return;
    }
    for (int round = 0; round < firstRuns; ++round) {
        if (!flatRun(mesh, state, hard, 0, round, -1, -1).again) {
            break;
        }
    }

    const int budget = state.refineVertexBudget < 1
                           ? truncToInt(static_cast<float>(PrepState::workResolution) * budgetFactor)
                           : state.refineVertexBudget;
    EdgeSplitter splitter(passes == -1 ? defaultPasses : passes, budget, state.freezeBorders);
    // Finding the symmetry splits mirror edges together and fixes the
    // diagonal of triangles with two split edges.
    splitter.setMirrorSplit(state.symmetry.detect);
    splitter.setShortestDiagonal(!state.symmetry.detect);
    const float lr = targetLength(mesh, state, area);
    if (lr > 0.0f) {
        splitter.begin();
        const BitArray* locked = (PrepState::useGroupBorders || state.hasEdgeAttribute) ? &hard : nullptr;
        while (true) {
            const int p = splitter.passCount();
            if (p >= 2 && !splitter.needsAnotherPass()) {
                break;
            }
            if (p > 3) {
                break;
            }
            splitter.pass(mesh, (p == 0 || p == 3) ? lr + lr : lr, locked);
            if (splitter.passInserted() > 0) {
                // The hard edges follow the new edges, lines over up to two edges.
                hard = state.lines.collectEdges(mesh, PrepState::useGroupBorders, state.hasEdgeAttribute, true,
                                                NetLines::Restore::OverTwoEdges);
            }
            const int q = splitter.passCount();
            bool check = true;
            if (q > 1) {
                if (!splitter.needsAnotherPass() || q > 3) {
                    break;
                }
                check = splitter.passFirstNew() < splitter.passSharpEnd();
            }
            if (check) {
                (void)flatRun(mesh, state, hard, 1, p, splitter.passFirstNew(), splitter.passSharpEnd());
            }
        }
    }
    state.refineInsertedVertexCount = splitter.insertedTotal();
    (void)flatRun(mesh, state, hard, 2, 0, -1, -1);
}

} // namespace areko::remesher::prep
