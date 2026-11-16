#pragma once

#include "trace.h"
#include "trace_frame.h"

namespace areko::remesher::extract {

// Finding a node again at a place, and putting a node onto a mesh edge.
class NodeSearch {
public:
    // A node within the squared distance `toleranceSq` of `position`: first in
    // the list of triangle t (newest first); with `extended` then in the
    // triangle across edge `side` of t (side >= 0; across a cut edge with the
    // tolerance for cut edges of the trace), then in every triangle at the
    // three points of t that was not searched yet. -1 when there is none.
    static int find(const TraceFrame& frame, const Trace& trace, float toleranceSq, int t, int side,
                    mesh::Vec3f position, bool extended);

    // Makes `node` a node on an edge of its triangle when it lies closer to
    // that edge than a hundred thousandth of the doubled triangle area
    // (squared distance); the nearest edge wins, the earlier one on a tie.
    static void snapToEdge(const TraceFrame& frame, TraceNode& node);

    // Squared distance from p to the segment a-b: the end a when p lies
    // behind it, else the foot on the line when p lies before b (or the test
    // is not decided), else the end b.
    [[nodiscard]] static float segmentDistanceSq(mesh::Vec3f p, mesh::Vec3f a, mesh::Vec3f b);
};

} // namespace areko::remesher::extract
