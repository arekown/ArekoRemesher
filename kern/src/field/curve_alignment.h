#pragma once

#include "cross_field.h"
#include "cut_graph.h"

#include "mesh/mesh.h"
#include "mesh/vertex_fan.h"

#include <cstddef>
#include <vector>

namespace areko::remesher::field {

// What a walk around a point of a curve reads from the combed cross field,
// and the walk itself: the fan of the point, the direction label of a curve
// edge seen from one of its triangles, how a cut edge turns the label, and
// whether a sector of the fan bends back. The mesh, the field and the cut
// graph must outlive it.
class CurveLabels {
public:
    // Where the walk of a caller departs from the plain one.
    struct WalkRules {
        // A fan that is empty because its start edge is a border edge at the
        // far end of the fan is taken from the other border edge of the point.
        bool startAtOtherBorderEdge = false;
        // The walk ends in front of the curve edge its label was last taken at.
        bool endBeforeReference = false;
        // A curve check cut between two curve edges separates them already:
        // the second one is not compared.
        bool checkCutsSeparate = false;
    };

    CurveLabels(const mesh::Mesh& mesh, const CrossField& field, const CutGraph& graph);

    // The end points of all curve edges.
    [[nodiscard]] mesh::BitArray curvePoints() const;
    // Label of the curve direction from point `from` to point `to`, seen in
    // the triangle of `corner`: 0 along the first combed direction, 2 against
    // it, 1 along the second, 3 against it. A product that is not a number
    // counts as against.
    [[nodiscard]] int label(int corner, int from, int to) const;

    // Walks the fan around point v from its first curve edge (in the edge
    // order of v) and carries the label of every curve edge through the
    // triangles to the next one: the label is taken behind a curve edge,
    // looking away from v; a cut edge on the way turns it. A curve edge that
    // is reached with the label it has itself, seen from the triangle in
    // front of it, runs the same way as the curve edge before it:
    // `sameLabel(fan, from, to)` is called with the fan positions of the two,
    // unless the sector between them bends back. Behind the last edge of the
    // fan no label is taken. Nothing happens at a point without curve edge or
    // with an empty fan.
    template <class SameLabel>
    void walkAround(int v, const WalkRules& rules, SameLabel sameLabel) const;

private:
    // First curve edge in the edge order of point v, or -1.
    [[nodiscard]] int firstCurveEdgeAt(int v) const;
    // Fan around v from `startEdge`, walked to its end however many corners
    // the point has.
    [[nodiscard]] mesh::VertexFan fanFrom(int v, int startEdge) const
    {
        return mesh::vertexFan(adjacency_, v, startEdge, false);
    }
    // The label behind edge `edge`, reached from the triangle of corner
    // `before` (the corner of that triangle the edge starts or ends at): a
    // cut edge turns the label by the quarter turns left there, any other
    // edge leaves it.
    [[nodiscard]] int behind(int label, int edge, int before) const;
    // Whether the sector from fan corner i up to fan edge j is wider than a
    // straight angle and one of its inner edges points away from edge j.
    [[nodiscard]] bool bendsBack(int v, const mesh::VertexFan& fan, std::size_t i, std::size_t j) const;

    const mesh::Adjacency& adjacency_;
    const std::vector<mesh::Vec3f>& points_;
    const CrossField& field_;
    const CutGraph& graph_;
    std::vector<float> cornerAngles_;
};

template <class SameLabel>
void CurveLabels::walkAround(int v, const WalkRules& rules, SameLabel sameLabel) const
{
    const int startEdge = firstCurveEdgeAt(v);
    if (startEdge == -1) {
        return;
    }
    mesh::VertexFan fan = fanFrom(v, startEdge);
    if (rules.startAtOtherBorderEdge && fan.corners.empty() && adjacency_.isBoundaryEdge(startEdge)) {
        if (const int other = adjacency_.firstBoundaryEdgeAt(v, startEdge); other != -1) {
            fan = fanFrom(v, other);
        }
    }
    if (fan.corners.empty()) {
        return;
    }
    const std::vector<int>& corners = fan.corners;
    const std::vector<int>& edges = fan.edges;
    int current = label(corners[0], adjacency_.vertex(corners[0]), adjacency_.vertex(adjacency_.next(corners[0])));
    int reference = edges[0];
    bool crossedCheckCut = false;
    std::size_t from = 0;
    for (std::size_t i = 1; i < edges.size(); ++i) {
        const int e = edges[i];
        if (rules.endBeforeReference && e == reference) {
            return;
        }
        const int before = adjacency_.prev(corners[i - 1]);
        if (rules.checkCutsSeparate && graph_.curveCheckCutBits.testBounded(static_cast<std::size_t>(e))) {
            crossedCheckCut = true;
        }
        if (graph_.curveEdgeKind[static_cast<std::size_t>(e)] == curve_edge::none) {
            current = behind(current, e, before);
            continue;
        }
        if (!crossedCheckCut) {
            const int arriving = label(before, adjacency_.vertex(corners[i - 1]), adjacency_.vertex(before));
            if (current == arriving && !bendsBack(v, fan, from, i)) {
                sameLabel(fan, from, i);
            }
        }
        if (i != edges.size() - 1) {
            current = label(corners[i], adjacency_.vertex(corners[i]),
                            adjacency_.vertex(adjacency_.next(corners[i])));
            reference = e;
            crossedCheckCut = false;
            from = i;
        }
    }
}

// How the curve edges lie in the combed cross field.
//  - every half edge of a curve edge follows the first or the second combed
//    direction of its triangle, whichever it runs along more closely,
//  - around every point of a curve the direction label of one curve edge is
//    carried through the triangles to the next curve edge (four labels: along
//    or against the first or the second direction; a cut edge on the way
//    turns the label by its residual jump); two neighbouring curve edges that
//    end up with the same label cannot both become edge loops of a quad
//    mesh, unless the sector between them is wider than a straight angle and
//    bends back,
//  - for such a pair one edge between them becomes a curve check cut and the
//    point a curve check source of the tree, so that the next cut graph
//    separates them.
class CurveAlignment {
public:
    // Whether the half edge that starts at `corner` runs along the first
    // combed direction of its triangle rather than the second. A product
    // that is not a number counts as the second.
    [[nodiscard]] static bool followsFirst(const mesh::Mesh& mesh, const CrossField& field, const CutGraph& graph,
                                           int corner);

    // Writes curveFollowsU of `out` (one bit per corner).
    static void orient(const mesh::Mesh& mesh, const CrossField& field, CutGraph& out);

    // Appends to curveCheckCuts, curveCheckCutBits and curveCheckSources of `out`.
    static void check(const mesh::Mesh& mesh, const CrossField& field, CutGraph& out);
};

} // namespace areko::remesher::field
