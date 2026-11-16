#include "trace_frame.h"

#include "mesh/vertex_fan.h"

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

TraceFrame::TraceFrame(const parts::Part& part, Trace& trace)
    : mesh_(part.mesh.adjacency()), copy_(part.cutGraph->copy), uvMesh_(part.cutGraph->copy.adjacency()),
      uv_(part.solution->uv), points_(part.mesh.points()), doubleAreas_(part.mesh.triangleGeometry().doubleAreas()),
      residualJump_(part.cutGraph->residualJump), cutEdges_(part.cutGraph->cutEdgeBits),
      frozen_(part.frozenBorderMarks()), triangles_(part.mesh.faceCount()), flipState_(trace.flipState)
{
    flipState_.resize(static_cast<std::size_t>(triangles_), 0);
}

bool TraceFrame::borderPoint(int point) const
{
    for (const int e : mesh_.pointEdges(point)) {
        if (mesh_.isBoundaryEdge(e)) {
            return true;
        }
    }
    return false;
}

bool TraceFrame::seamless(int e) const
{
    return copy_.copyEdges(mesh_, e).size() == 1;
}

bool TraceFrame::flipped(int t)
{
    std::uint8_t& state = flipState_[static_cast<std::size_t>(t)];
    if (state == 0) {
        const std::array<float, 2>& p0 = cornerUv(3 * t);
        const std::array<float, 2>& p1 = cornerUv(3 * t + 1);
        const std::array<float, 2>& p2 = cornerUv(3 * t + 2);
        const float area = (p2[1] - p0[1]) * (p1[0] - p0[0]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
        state = area < 0.0f ? 1 : 2;
    }
    return state == 1;
}

int TraceFrame::cornerAt(int t, int point) const
{
    for (int c = 3 * t; c < 3 * t + 3; ++c) {
        if (mesh_.vertex(c) == point) {
            return c;
        }
    }
    return -1;
}

bool TraceFrame::aroundPoint(int from, int to, int& s, bool forward) const
{
    const VertexFan fan = vertexFanFrom(mesh_, from, forward);
    // Step j crosses edge j + 1 from corner j to corner j + 1. A fan of more
    // than 2000 corners has no edges and takes no step.
    for (std::size_t j = 0; j + 1 < fan.corners.size() && j + 1 < fan.edges.size(); ++j) {
        const int g = fan.edges[j + 1];
        int r = residualJump(g);
        if (fan.corners[j] / 3 != mesh_.edgeCorners(g)[0] / 3) {
            r = -r;
        }
        s = mod4(s - r);
        if (fan.corners[j + 1] == to) {
            return true;
        }
    }
    return false;
}

std::optional<int> TraceFrame::directionSlot(int t1, const TraceNode& node, int d, int t2)
{
    if (t1 == t2) {
        return d;
    }
    int s = d;
    if (uvMesh_.commonEdge(t1, t2) == -1) {
        const int e = mesh_.commonEdge(t1, t2);
        if (e >= 0) {
            int r = residualJump(e);
            if (t2 != mesh_.edgeCorners(e)[0] / 3) {
                r = -r;
            }
            s = mod4(d - r);
        } else {
            // Only a point in common: sum the turns around it.
            int q = -1;
            if (node.kind == node_kind::point) {
                q = node.reference;
            } else {
                for (int c = 3 * t1; c < 3 * t1 + 3 && q < 0; ++c) {
                    if (cornerAt(t2, mesh_.vertex(c)) >= 0) {
                        q = mesh_.vertex(c);
                    }
                }
            }
            if (q < 0) {
                return std::nullopt;
            }
            const int i1 = cornerAt(t1, q);
            const int i2 = cornerAt(t2, q);
            if (i1 < 0 || i2 < 0) {
                return std::nullopt;
            }
            if (uvMesh_.vertex(i1) == uvMesh_.vertex(i2)) {
                return d;
            }
            int forward = d;
            if (aroundPoint(i2, i1, forward, true)) {
                return forward;
            }
            // The fan the other way round starts again from d.
            int backward = d;
            if (aroundPoint(i2, i1, backward, false)) {
                return backward;
            }
            return std::nullopt;
        }
    }
    if (flipped(t1) != flipped(t2)) {
        s = (s + 2) % 4;
    }
    return s;
}

} // namespace areko::remesher::extract
