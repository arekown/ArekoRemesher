#pragma once

#include "segment_tracer.h"
#include "trace.h"
#include "trace_frame.h"

#include <array>

namespace areko::remesher::extract {

// Walks around one face of the grid: from a node along one direction to the
// next node, there always turning the same way (a quarter turn to the left,
// to the right in a turned-over triangle), until the walk comes back to its
// start. A walk that closes is handed to the face closer.
class FaceWalker {
public:
    FaceWalker(TraceFrame& frame, Trace& trace, SegmentTracer& segments, FaceCloser& closer)
        : frame_(frame), trace_(trace), segments_(segments), closer_(closer)
    {
    }

    // Walk from node `start` in `direction` in triangle `triangle` with grid
    // coordinates `uv`; true when the face closed.
    bool walk(std::array<float, 2> uv, int start, int direction, int triangle);

    // Whether a polygon of the same length as `loop` holds the numbers of
    // its nodes in the same order from the same place; only loops of seed
    // nodes count.
    [[nodiscard]] static bool faceExists(const Trace& trace, const FaceLoop& loop);

    // Whether the loop, just grown by its last node, is to be given up;
    // `unchanged` counts the steps that left the coordinates as they were.
    [[nodiscard]] static bool tooLong(const FaceLoop& loop, int unchanged);

private:
    // Gives every node of the loop its output number and marks the
    // direction it was left in as done.
    void settle(const FaceLoop& loop);

    TraceFrame& frame_;
    Trace& trace_;
    SegmentTracer& segments_;
    FaceCloser& closer_;
};

} // namespace areko::remesher::extract
