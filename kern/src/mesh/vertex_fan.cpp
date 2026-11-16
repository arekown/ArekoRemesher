#include "vertex_fan.h"

namespace areko::remesher::mesh {

namespace {

constexpr std::size_t maxCorners = 2000;

// Walk around v from corner `first`. Forward (`atPoint` true) each step
// crosses the edge of the previous corner, backward the edge of the corner
// itself. `atPoint` ends as the direction of the last step taken.
VertexFan ring(const Adjacency& adjacency, int v, int first, bool& atPoint, bool limited = true)
{
    VertexFan fan{{first}, {adjacency.cornerEdge(atPoint ? first : adjacency.prev(first))}};
    int c = first;
    while (true) {
        if (atPoint) {
            c = adjacency.prev(c);
        }
        const int e = adjacency.cornerEdge(c);
        if (e < 0) {
            break;
        }
        int other = adjacency.otherCorner(e, c);
        if (other < 0) {
            fan.edges.push_back(e);
            break;
        }
        atPoint = adjacency.vertex(other) == v;
        if (!atPoint) {
            other = adjacency.next(other);
        }
        if (other == first) {
            fan.edges.push_back(e);
            break;
        }
        if (limited && fan.corners.size() > maxCorners) {
            const auto all = adjacency.pointCorners(v);
            fan.corners.insert(fan.corners.end(), all.begin(), all.end());
            fan.edges.clear();
            break;
        }
        fan.edges.push_back(e);
        fan.corners.push_back(other);
        c = other;
    }
    return fan;
}

VertexFan forwardRing(const Adjacency& adjacency, int v, int first, bool limited = true)
{
    bool atPoint = true;
    return ring(adjacency, v, first, atPoint, limited);
}

} // namespace

VertexFan vertexFan(const Adjacency& adjacency, int v, int startEdge, bool limited)
{
    const auto corners = adjacency.edgeCorners(startEdge);
    if (corners.size() == 1) {
        return adjacency.vertex(corners[0]) == v ? forwardRing(adjacency, v, corners[0], limited) : VertexFan{};
    }
    if (corners.size() == 2) {
        if (adjacency.vertex(corners[0]) == v) {
            return forwardRing(adjacency, v, corners[0], limited);
        }
        if (adjacency.vertex(corners[1]) == v) {
            return forwardRing(adjacency, v, corners[1], limited);
        }
    }
    return {};
}

VertexFan vertexFanFrom(const Adjacency& adjacency, int corner, bool forward)
{
    bool atPoint = forward;
    return ring(adjacency, adjacency.vertex(corner), corner, atPoint);
}

VertexFan orderedVertexFan(const Adjacency& adjacency, int v)
{
    const auto corners = adjacency.pointCorners(v);
    if (corners.empty()) {
        return {};
    }
    const int first = corners[0];
    if (corners.size() == adjacency.pointEdges(v).size()) {
        VertexFan fan = forwardRing(adjacency, v, first);
        if (fan.corners.size() == corners.size()) {
            return fan;
        }
    }
    bool atPoint = false;
    const VertexFan back = ring(adjacency, v, first, atPoint);
    atPoint = !atPoint;
    return ring(adjacency, v, back.corners.back(), atPoint);
}

std::vector<int> vertexFanEdges(const Adjacency& adjacency, int v, int startEdge)
{
    return vertexFan(adjacency, v, startEdge).edges;
}

} // namespace areko::remesher::mesh
