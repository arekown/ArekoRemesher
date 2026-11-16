#include "isoline_tracer.h"

#include "face_walker.h"
#include "segment_tracer.h"
#include "trace_frame.h"

#include "param/solution_seeds.h"

#include <utility>
#include <vector>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

// Seed nodes for `seeds`, one per point (-1 for a point without corner):
// on the point, numbered at once, not on the stack.
std::vector<int> seedNodes(const TraceFrame& frame, TraceNodes& nodes, const std::vector<int>& seeds)
{
    std::vector<int> out;
    out.reserve(seeds.size());
    for (const int p : seeds) {
        const auto corners = frame.mesh().pointCorners(p);
        if (corners.empty()) {
            out.push_back(-1);
            continue;
        }
        const int k = nodes.create(corners[0] / 3, frame.position(p), {0.0f, 0.0f}, true);
        nodes.dropTop();
        TraceNode& node = nodes[k];
        node.triangle = -1;
        node.kind = node_kind::point;
        node.reference = p;
        node.flags = node_flag::seed;
        if (frame.borderPoint(p)) {
            node.flags |= node_flag::border;
        }
        node.id = nodes.takeId();
        out.push_back(k);
    }
    return out;
}

} // namespace

void IsolineTracer::trace(parts::Part& part, FaceCloser& closer, const run::Progress& progress)
{
    part.trace.reset();
    if (!part.cutGraph || !part.solution) {
        return;
    }
    const field::CutGraph& graph = *part.cutGraph;
    const param::Solution& solution = *part.solution;
    if (graph.residualJump.empty() || (solution.singularities.points.empty() && solution.seeds.empty())) {
        return;
    }
    const int triangles = part.mesh.faceCount();
    Trace result(triangles);
    result.segmentVisits.assign(static_cast<std::size_t>(triangles), 0);
    if (!part.curves.empty()) {
        result.curves = CurveMarks::mark(part);
    }
    TraceFrame frame(part, result);
    SegmentTracer segments(frame, result);
    FaceWalker walker(frame, result, segments, closer);
    TraceNodes& nodes = result.nodes;

    // Walks from the start points; once more with the border points added
    // when that closes no face.
    const bool fromSingular = !solution.singularities.points.empty();
    std::vector<int> seeds = fromSingular ? solution.singularities.points : solution.seeds;
    bool reseeded = false;
    for (;;) {
        const TraceNodes::Mark firstPass = nodes.mark();
        const std::vector<int> starts = seedNodes(frame, nodes, seeds);
        const float diagonal = mesh::diagonal(part.mesh.bounds());
        result.nodeTolerance = diagonal * 1e-7f;
        result.nodeToleranceSq = result.nodeTolerance * result.nodeTolerance;
        result.markedTolerance = diagonal * 1e-6f;
        result.markedToleranceSq = result.markedTolerance * result.markedTolerance;
        segments.setStackPhase(false);
        for (std::size_t i = 0; i < starts.size(); ++i) {
            const int k = starts[i];
            if (k < 0) {
                continue;
            }
            if (i % 50 == 0) {
                progress.check();
            }
            for (const int q : graph.copy.pointCopies(frame.mesh(), seeds[i])) {
                const std::array<float, 2> uv = solution.uv[static_cast<std::size_t>(q)];
                if (!isWholeNumber(uv[0]) || !isWholeNumber(uv[1])) {
                    continue;
                }
                for (const int c : frame.uvMesh().pointCorners(q)) {
                    walker.walk(uv, k, 0, c / 3);
                    walker.walk(uv, k, 1, c / 3);
                }
            }
            if (!fromSingular) {
                nodes[k].flags |= node_flag::directions;
            }
        }
        if (!reseeded && result.polygons.rows() == 0) {
            const std::size_t before = seeds.size();
            param::SolutionSeeds::appendBorderPoints(graph.copy, part.mesh.pointCount(), solution.uv, seeds);
            if (seeds.size() > before) {
                reseeded = true;
                // The nodes of the first pass go before it starts again,
                // and so do its counts.
                nodes.discardSince(firstPass);
                result.counts = {};
                result.segmentVisits.assign(static_cast<std::size_t>(triangles), 0);
                continue;
            }
        }
        break;
    }

    // Walks from every node made on the way, newest first.
    progress.check();
    segments.setStackPhase(true);
    while (!nodes.stackEmpty()) {
        const int k = nodes.pop();
        for (int d = 0; d < 4; ++d) {
            if (nodes[k].done(d)) {
                continue;
            }
            const bool closed = walker.walk(nodes[k].uv, k, d, nodes[k].triangle);
            // A stop request is looked for at every 500th face.
            if (closed && result.polygons.rows() % 500 == 0) {
                progress.check();
            }
        }
    }
    segments.setStackPhase(false);

    // The walks put back in that phase.
    const std::size_t deferred = result.deferredNodes.size();
    for (std::size_t j = 0; j < deferred; ++j) {
        const int k = result.deferredNodes[j];
        const int d = result.deferredDirections[j];
        if (!nodes[k].done(d)) {
            walker.walk(nodes[k].uv, k, d, nodes[k].triangle);
        }
    }
    part.trace = std::move(result);
}

} // namespace areko::remesher::extract
