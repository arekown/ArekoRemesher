#include "node_search.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

int searchIn(const Trace& trace, int t, float toleranceSq, Vec3f position)
{
    return trace.nodes.findIn(t, [&](int n) {
        const Vec3f d = trace.nodes[n].position - position;
        return dot(d, d) < toleranceSq;
    });
}

} // namespace

int NodeSearch::find(const TraceFrame& frame, const Trace& trace, float toleranceSq, int t, int side,
                     Vec3f position, bool extended)
{
    int found = searchIn(trace, t, toleranceSq, position);
    if (found >= 0 || !extended) {
        return found;
    }
    const Adjacency& mesh = frame.mesh();
    std::vector<int> searched{t};
    if (side >= 0) {
        const int c = 3 * t + side;
        const int e = mesh.cornerEdge(c);
        const int other = e >= 0 ? mesh.otherCorner(e, c) : -1;
        if (other >= 0) {
            const int across = other / 3;
            const float tolerance = frame.cutEdge(e) ? trace.markedToleranceSq : toleranceSq;
            found = searchIn(trace, across, tolerance, position);
            if (found >= 0) {
                return found;
            }
            searched.push_back(across);
        }
    }
    for (int c = 3 * t; c < 3 * t + 3; ++c) {
        for (const int corner : mesh.pointCorners(mesh.vertex(c))) {
            const int n = corner / 3;
            if (std::find(searched.begin(), searched.end(), n) != searched.end()) {
                continue;
            }
            found = searchIn(trace, n, toleranceSq, position);
            if (found >= 0) {
                return found;
            }
            searched.push_back(n);
        }
    }
    return -1;
}

void NodeSearch::snapToEdge(const TraceFrame& frame, TraceNode& node)
{
    const int t = node.triangle;
    float limit = frame.doubleArea(t) / 100000.0f;
    for (int c = 3 * t; c < 3 * t + 3; ++c) {
        const Vec3f a = frame.position(frame.cornerPoint(c));
        const Vec3f b = frame.position(frame.cornerPoint(frame.mesh().next(c)));
        const float d2 = segmentDistanceSq(node.position, a, b);
        const int e = frame.mesh().cornerEdge(c);
        if (d2 < limit && e != -1) {
            node.kind = node_kind::edge;
            node.reference = e;
            limit = d2;
        }
    }
}

float NodeSearch::segmentDistanceSq(Vec3f p, Vec3f a, Vec3f b)
{
    return closestOnSegment(p, a, b).distanceSq;
}

} // namespace areko::remesher::extract
