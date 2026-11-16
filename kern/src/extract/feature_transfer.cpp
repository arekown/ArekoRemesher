#include "feature_transfer.h"

#include "curve_cut_transfer.h"
#include "group_transfer.h"
#include "hard_edge_transfer.h"

#include "parts/part.h"
#include "util/introsort.h"

#include <numeric>
#include <vector>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

// Output points that come from a point or an edge of a crease chain.
void markOnCurve(const CurveMarks& marks, OutputPoints& points)
{
    for (std::size_t i = 0; i < points.source.size(); ++i) {
        const PointSource& s = points.source[i];
        if (s.id < 0) {
            continue;
        }
        const auto id = static_cast<std::size_t>(s.id);
        if ((s.type == 1 && marks.points.testBounded(id)) || (s.type == 2 && marks.edges.testBounded(id))) {
            points.onCurve.setBounded(i);
        }
    }
}

// Every marked point gives a third of its weight to each neighbour that is
// not marked (at most 1 there) and keeps 0.
void spreadWeights(const Adjacency& adjacency, OutputPoints& points)
{
    points.onCurve.forEachSet([&](std::size_t i) {
        const float w = points.weight[i];
        points.weight[i] = 0.0f;
        if (!(0.0f < w)) {
            return;
        }
        for (const int j : adjacency.pointNeighbours(static_cast<int>(i))) {
            if (points.onCurve.testBounded(static_cast<std::size_t>(j))) {
                continue;
            }
            float& to = points.weight[static_cast<std::size_t>(j)];
            const float x = w / 3.0f + to;
            to = 1.0f <= x ? 1.0f : x;
        }
    });
}

// One polyline per crease chain: the marked output points with an edge, at
// their place along the chain, sorted by it; a closed chain repeats its first
// entry at the end.
std::vector<CurvePolyline> buildPolylines(const Mesh& triangles, const CurveMarks& marks, const RawMesh& raw)
{
    std::vector<CurvePolyline> lines(static_cast<std::size_t>(marks.chainEdges.rows()));
    const Adjacency& output = raw.mesh.adjacency();
    const Adjacency& input = triangles.adjacency();
    const auto& inputPoints = triangles.points();
    const auto& outputPoints = raw.mesh.points();
    auto add = [&](int chain, std::size_t i, float key) {
        lines[static_cast<std::size_t>(chain)].points.push_back(static_cast<int>(i));
        lines[static_cast<std::size_t>(chain)].keys.push_back(key);
    };
    raw.points.onCurve.forEachSet([&](std::size_t i) {
        if (output.pointEdges(static_cast<int>(i)).empty()) {
            return;
        }
        const PointSource& s = raw.points.source[i];
        if (s.type == 1) {
            for (const int e : input.pointEdges(s.id)) {
                const int chain = marks.edgeChain[static_cast<std::size_t>(e)];
                if (chain == -1) {
                    continue;
                }
                const bool along = marks.alongFirstEnd.testBounded(static_cast<std::size_t>(e));
                const bool first = input.edgeVertices(e)[0] == s.id;
                const float offset = first == along ? 0.0f : 1.0f;
                add(chain, i, static_cast<float>(marks.edgePlace[static_cast<std::size_t>(e)]) + offset);
            }
        } else if (s.type == 2) {
            const int e = s.id;
            const int chain = marks.edgeChain[static_cast<std::size_t>(e)];
            if (chain == -1) {
                return;
            }
            const auto ends = input.edgeVertices(e);
            const Vec3f a = inputPoints[static_cast<std::size_t>(ends[0])];
            const float d1 = distance(a, outputPoints[i]);
            const float d2 = distance(inputPoints[static_cast<std::size_t>(ends[1])], a);
            const float r = d1 / d2;
            float t = 1.0f < r ? 1.0f : r;
            if (!marks.alongFirstEnd.testBounded(static_cast<std::size_t>(e))) {
                t = 1.0f - t;
            }
            add(chain, i, t + static_cast<float>(marks.edgePlace[static_cast<std::size_t>(e)]));
        }
    });
    std::vector<int> order;
    for (std::size_t c = 0; c < lines.size(); ++c) {
        CurvePolyline& line = lines[c];
        order.resize(line.points.size());
        std::iota(order.begin(), order.end(), 0);
        util::Introsort::sort(order.begin(), order.end(), [&](int x, int y) {
            return line.keys[static_cast<std::size_t>(x)] < line.keys[static_cast<std::size_t>(y)];
        });
        CurvePolyline sorted;
        for (const int k : order) {
            sorted.points.push_back(line.points[static_cast<std::size_t>(k)]);
            sorted.keys.push_back(line.keys[static_cast<std::size_t>(k)]);
        }
        const auto chain = marks.chainPoints.row(static_cast<int>(c));
        if (!chain.empty() && chain.front() == chain.back() && !sorted.points.empty()) {
            sorted.points.push_back(sorted.points.front());
            sorted.keys.push_back(sorted.keys.front());
        }
        line = std::move(sorted);
    }
    return lines;
}

} // namespace

void FeatureTransfer::run(parts::Part& part, bool cutAlongCurves)
{
    if (!part.raw || !part.trace) {
        return;
    }
    RawMesh& raw = *part.raw;
    const CurveMarks& marks = part.trace->curves;
    raw.points.onCurve = BitArray(static_cast<std::size_t>(raw.mesh.pointCount()));
    if (marks.edgeCount == 0) {
        GroupTransfer::uniform(part.mesh, raw);
        return;
    }
    markOnCurve(marks, raw.points);
    spreadWeights(raw.mesh.adjacency(), raw.points);
    const std::vector<CurvePolyline> lines = buildPolylines(part.mesh, marks, raw);
    if (cutAlongCurves) {
        CurveCutTransfer::run(part.mesh, marks, lines, raw);
    } else {
        const Adjacency& adjacency = raw.mesh.adjacency();
        raw.curveEdges = BitArray(static_cast<std::size_t>(adjacency.edgeCount()));
        for (const CurvePolyline& line : lines) {
            CurveCutTransfer::markEdges(adjacency, line.points, raw.curveEdges);
        }
    }
    HardEdgeTransfer::run(part.mesh, part.curves, raw);
    GroupTransfer::run(part.mesh, raw);
}

} // namespace areko::remesher::extract
