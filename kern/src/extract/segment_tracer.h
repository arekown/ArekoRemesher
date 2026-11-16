#pragma once

#include "trace.h"
#include "trace_frame.h"

#include <array>
#include <optional>
#include <vector>

namespace areko::remesher::extract {

// Follows one grid line from a node, triangle by triangle, until the moving
// coordinate reaches the next whole value: the fixed coordinate stays on its
// whole value, the moving one runs in direction d (0: u grows, 1: v grows,
// 2: u falls, 3: v falls). Across a seam of the cut copy the coordinates and
// the direction are carried over with the quarter turns of the seam. At the
// end a node is found again or made, and both nodes remember the segment.
class SegmentTracer {
public:
    SegmentTracer(TraceFrame& frame, Trace& trace) : frame_(frame), trace_(trace) {}

    // code: 0 done, 1 ran into the border, 2 left the band of +-1 around the
    // start value, 3 looped or crossed too many edges, 4 no way out of a
    // triangle, 5 invalid, 6 put back (stack phase, the segment came back to
    // its own grid line). With code 0: the node reached (-1 only when a seed
    // node has no corner in its first triangle), the direction of arrival,
    // and whether the segment came back to its own grid line.
    struct Result {
        int code = 0;
        int node = -1;
        int arrival = 0;
        bool degenerate = false;
    };

    // `triangle` and `uv` are where the segment starts and are moved to
    // where it ends. `firstStep`: the first segment of a walk from a seed.
    Result run(int start, int direction, int& triangle, std::array<float, 2>& uv, bool firstStep);

    // While the stack phase runs, a segment that comes back to its own grid
    // line is put back (code 6).
    void setStackPhase(bool on) noexcept { stackPhase_ = on; }
    // Collect the kept border points the segment passes as nodes (cleared
    // before every segment by the caller).
    void setCollecting(bool on) noexcept { collecting_ = on; }
    [[nodiscard]] bool collecting() const noexcept { return collecting_; }
    [[nodiscard]] std::vector<int>& collected() noexcept { return collected_; }

    // Position on the mesh of grid coordinates `uv` in triangle t:
    // barycentric in (u, v); when that is not a number, the point of the
    // corner the line leaves through (corner >= 0), the centre of a triangle
    // whose three corners have the same (u, v), else a point along the
    // longest side in (u, v). False when the replacement is not a number
    // either.
    static bool place(const TraceFrame& frame, int t, int corner, const std::array<float, 2>& uv,
                      mesh::Vec3f& out);

private:
    // The state of one segment while it runs from triangle to triangle.
    struct Walk;
    // The triangle the segment is in and the way out of it.
    struct Cell;
    // The whole value the moving coordinate reached.
    struct Target;

    Result follow(int start, int direction, int& triangle, std::array<float, 2>& uv, bool firstStep);
    // The corners of the current triangle in (u, v) and where the line
    // crosses its edges.
    [[nodiscard]] Cell enter(const Walk& walk) const;
    // The highest and the lowest crossing with an edge, with the fallbacks
    // when there is none, and which of the two the line leaves through.
    std::optional<Result> wayOut(Walk& walk, Cell& cell) const;
    // Whether the way out is a corner of the triangle.
    std::optional<Result> throughCorner(Walk& walk, Cell& cell) const;
    // Whether the moving coordinate reached a whole value on the way out.
    [[nodiscard]] static std::optional<Target> reached(const Walk& walk, const Cell& cell);
    // The end of the segment: its place, the direction of arrival, the node
    // found again or made, and the links of both nodes.
    Result finish(Walk& walk, const Cell& cell, const Target& target);
    // A new end node on `position` with its kind and border flag.
    int newNode(const Walk& walk, const Cell& cell, mesh::Vec3f position);
    // Both nodes remember the segment.
    void connect(const Walk& walk, int end, bool degenerate);
    // A kept border point of a frozen border the segment passes over a corner.
    void collectBorderPoint(const Walk& walk);
    // Into the next triangle, over the corner or the edge of the way out
    // (and on along the edge the line follows).
    std::optional<Result> crossOver(Walk& walk, Cell& cell) const;
    // Sets the border flag of a new node `node` that lies on a border edge of
    // its triangle along the fixed coordinate (axis `a` moving).
    void borderAlongLine(int node, int a);

    TraceFrame& frame_;
    Trace& trace_;
    bool stackPhase_ = false;
    bool collecting_ = false;
    std::vector<int> collected_;
};

} // namespace areko::remesher::extract
