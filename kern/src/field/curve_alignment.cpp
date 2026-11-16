#include "curve_alignment.h"

#include "mesh/measures.h"
#include "mesh/vertex_fan.h"

#include <algorithm>
#include <cmath>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

CurveLabels::CurveLabels(const Mesh& mesh, const CrossField& field, const CutGraph& graph)
    : adjacency_(mesh.adjacency()), points_(mesh.points()), field_(field), graph_(graph),
      cornerAngles_(cornerAngles(mesh))
{
}

BitArray CurveLabels::curvePoints() const
{
    BitArray curvePoints(static_cast<std::size_t>(adjacency_.pointCount()));
    for (const int e : graph_.curveEdges) {
        for (const int v : adjacency_.edgeVertices(e)) {
            curvePoints.set(static_cast<std::size_t>(v));
        }
    }
    return curvePoints;
}

int CurveLabels::firstCurveEdgeAt(int v) const
{
    for (const int e : adjacency_.pointEdges(v)) {
        if (graph_.curveEdgeKind[static_cast<std::size_t>(e)] != curve_edge::none) {
            return e;
        }
    }
    return -1;
}

int CurveLabels::label(int corner, int from, int to) const
{
    const std::size_t t = static_cast<std::size_t>(corner / 3);
    const bool alongU = graph_.curveFollowsU.testBounded(static_cast<std::size_t>(corner));
    const Vec3f direction = field_.inSpace(t, alongU ? graph_.combedU[t] : graph_.combedV[t]);
    const Vec3f d = points_[static_cast<std::size_t>(to)] - points_[static_cast<std::size_t>(from)];
    const bool against = !(dot(direction, d) > 0.0f);
    return alongU ? (against ? 2 : 0) : (against ? 3 : 1);
}

int CurveLabels::behind(int label, int edge, int before) const
{
    if (!graph_.cutEdgeBits.testBounded(static_cast<std::size_t>(edge))) {
        return label;
    }
    int jump = graph_.residualJump[static_cast<std::size_t>(edge)];
    if (before != adjacency_.edgeCorners(edge)[0]) {
        jump = -jump;
    }
    return ((label - jump) + 400) % 4;
}

bool CurveLabels::bendsBack(int v, const VertexFan& fan, std::size_t i, std::size_t j) const
{
    if (!(i + 1 < j)) {
        return false;
    }
    const Vec3f centre = points_[static_cast<std::size_t>(v)];
    const Vec3f last = points_[static_cast<std::size_t>(adjacency_.otherEnd(fan.edges[j], v))] - centre;
    float angle = cornerAngles_[static_cast<std::size_t>(fan.corners[i])];
    bool away = false;
    for (std::size_t k = i + 1; k < j; ++k) {
        angle = angle + cornerAngles_[static_cast<std::size_t>(fan.corners[k])];
        const Vec3f inner = points_[static_cast<std::size_t>(adjacency_.otherEnd(fan.edges[k], v))] - centre;
        if (dot(inner, last) < 0.0f) {
            away = true;
        }
    }
    return away && angle > piF;
}

namespace {

// The working state of one check around the points of the curves.
class Check {
public:
    Check(const Mesh& mesh, const CrossField& field, CutGraph& graph)
        : adjacency_(mesh.adjacency()), labels_(mesh, field, graph), graph_(graph)
    {
    }

    void run();

private:
    void separate(int v, const VertexFan& fan, std::size_t i, std::size_t j);

    const Adjacency& adjacency_;
    CurveLabels labels_;
    CutGraph& graph_;
};

void Check::separate(int v, const VertexFan& fan, std::size_t i, std::size_t j)
{
    const Adjacency& copy = graph_.copy.adjacency();
    for (std::size_t k = i + 1; k < j; ++k) {
        const int f = fan.edges[k];
        if (graph_.curveCheckCutBits.testBounded(static_cast<std::size_t>(f))) {
            continue;
        }
        // The other end must not lie on a cut or border already, and no
        // curve check cut may end there.
        const int w = adjacency_.otherEnd(f, v);
        const auto copyEdges = copy.pointEdges(w);
        if (std::any_of(copyEdges.begin(), copyEdges.end(), [&](int e) { return copy.isBoundaryEdge(e); })) {
            continue;
        }
        const auto edges = adjacency_.pointEdges(w);
        if (std::any_of(edges.begin(), edges.end(),
                        [&](int e) { return graph_.curveCheckCutBits.testBounded(static_cast<std::size_t>(e)); })) {
            continue;
        }
        if (graph_.curveCheckCutBits.empty()) {
            graph_.curveCheckCutBits = BitArray(static_cast<std::size_t>(adjacency_.edgeCount()));
        }
        graph_.curveCheckCutBits.set(static_cast<std::size_t>(f));
        graph_.curveCheckCuts.push_back(f);
        // Every point is a source once.
        std::vector<int>& sources = graph_.curveCheckSources;
        if (std::find(sources.begin(), sources.end(), v) == sources.end()) {
            sources.push_back(v);
        }
        return;
    }
}

void Check::run()
{
    // Two curve edges with the same label are separated when an edge lies between them.
    const CurveLabels::WalkRules rules{
        .startAtOtherBorderEdge = true, .endBeforeReference = true, .checkCutsSeparate = false};
    labels_.curvePoints().forEachSet([&](std::size_t point) {
        const int v = static_cast<int>(point);
        labels_.walkAround(v, rules, [&](const VertexFan& fan, std::size_t from, std::size_t to) {
            if (to - from > 1) {
                separate(v, fan, from, to);
            }
        });
    });
}

} // namespace

bool CurveAlignment::followsFirst(const Mesh& mesh, const CrossField& field, const CutGraph& graph, int corner)
{
    const Adjacency& adjacency = mesh.adjacency();
    const std::vector<Vec3f>& points = mesh.points();
    const std::size_t t = static_cast<std::size_t>(corner / 3);
    const Vec3f d = points[static_cast<std::size_t>(adjacency.vertex(adjacency.next(corner)))] -
                    points[static_cast<std::size_t>(adjacency.vertex(corner))];
    const float alongU = std::fabs(dot(field.inSpace(t, graph.combedU[t]), d));
    const float alongV = std::fabs(dot(field.inSpace(t, graph.combedV[t]), d));
    return alongU > alongV;
}

void CurveAlignment::orient(const Mesh& mesh, const CrossField& field, CutGraph& out)
{
    const Adjacency& adjacency = mesh.adjacency();
    out.curveFollowsU = BitArray(static_cast<std::size_t>(mesh.cornerCount()));
    for (const int e : out.curveEdges) {
        for (const int c : adjacency.edgeCorners(e)) {
            if (followsFirst(mesh, field, out, c)) {
                out.curveFollowsU.set(static_cast<std::size_t>(c));
            }
        }
    }
}

void CurveAlignment::check(const Mesh& mesh, const CrossField& field, CutGraph& out)
{
    Check(mesh, field, out).run();
}

} // namespace areko::remesher::field
