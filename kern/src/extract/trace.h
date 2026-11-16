#pragma once

#include "curve_marks.h"
#include "trace_nodes.h"

#include "mesh/csr.h"

#include <array>
#include <cstdint>
#include <vector>

namespace areko::remesher::extract {

// How the trace of a part went (without a first pass that was given up):
// segments by their outcome (0 reached the next grid line, also when a seed
// node has no corner in its first triangle, 1 ran into the border, 2 left
// the band of +-1, 3 looped or crossed too many edges, 4 found no way out
// of a triangle, 5 invalid, 6 put back), segments taken from the links of a
// node, and walks around faces.
struct TraceCounts {
    std::array<int, 7> segments{};
    int linkedSegments = 0;
    int walks = 0;
    int closedWalks = 0;
    // Walks given up at the limit of their length.
    int tooLong = 0;
};

// The result of the trace of one part. Later steps of the extraction read
// it and add to it.
struct Trace {
    explicit Trace(int triangleCount) : nodes(triangleCount) {}

    TraceNodes nodes;
    // Output polygons as node numbers, in the order the faces were closed.
    mesh::Csr polygons;
    CurveMarks curves;
    // Per triangle 0 unknown, 1 turned over in (u, v), 2 not turned over.
    std::vector<std::uint8_t> flipState;
    // Distance below which a node is found again (and its square), and the
    // same across a cut edge.
    float nodeTolerance = 0.0f;
    float nodeToleranceSq = 0.0f;
    float markedTolerance = 0.0f;
    float markedToleranceSq = 0.0f;
    // Walks put back during the stack phase (start node, start direction).
    std::vector<int> deferredNodes;
    std::vector<int> deferredDirections;
    TraceCounts counts;
    // Per triangle how many traced segments ran through it.
    std::vector<int> segmentVisits;
};

// A walk that came back to its start: its nodes, the direction followed from
// each, the triangle it was followed in, and per side whether its segment
// came back to the grid line it started on. The last side runs from the last
// node back to the first.
struct FaceLoop {
    std::vector<int> nodes;
    std::vector<int> directions;
    std::vector<int> triangles;
    std::vector<std::uint8_t> degenerate;
};

// Closes the faces the trace finds. The trace has already given every node
// of `loop` its output number and marked the directions it walked as done;
// the closer turns the loop into polygons of trace.polygons.
class FaceCloser {
public:
    FaceCloser() = default;
    FaceCloser(const FaceCloser&) = delete;
    FaceCloser& operator=(const FaceCloser&) = delete;
    virtual ~FaceCloser() = default;
    virtual void close(FaceLoop& loop, Trace& trace) = 0;
};

// Every loop becomes one polygon of the numbers of its nodes, as it is.
class LoopPolygons final : public FaceCloser {
public:
    void close(FaceLoop& loop, Trace& trace) override;
};

} // namespace areko::remesher::extract
