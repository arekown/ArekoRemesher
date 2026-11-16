#pragma once

#include "trace.h"
#include "trace_node.h"

#include "mesh/adjacency.h"
#include "mesh/attribute_map.h"
#include "parts/part.h"

#include <array>
#include <optional>
#include <span>
#include <vector>

namespace areko::remesher::extract {

// Remainder of x by 4, never negative.
[[nodiscard]] constexpr int mod4(int x) noexcept
{
    const int m = x % 4;
    return m < 0 ? m + 4 : m;
}

// What the trace of one part looks at: the triangle mesh of the part (its
// points, neighbourhood and doubled triangle areas), the cut copy with the
// grid coordinates of its points, the quarter turns left across every edge,
// the cut edges and the freeze marks; and whether a triangle is turned over
// in (u, v), worked out once per triangle.
class TraceFrame {
public:
    // The part must have a cut graph and a solution.
    TraceFrame(const parts::Part& part, Trace& trace);
    TraceFrame(const TraceFrame&) = delete;
    TraceFrame& operator=(const TraceFrame&) = delete;

    [[nodiscard]] const mesh::Adjacency& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const mesh::Adjacency& uvMesh() const noexcept { return uvMesh_; }
    [[nodiscard]] const field::CutCopy& copy() const noexcept { return copy_; }
    [[nodiscard]] int triangleCount() const noexcept { return triangles_; }

    // Grid coordinates and point of the mesh at corner c.
    [[nodiscard]] const std::array<float, 2>& cornerUv(int c) const
    {
        return uv_[static_cast<std::size_t>(uvMesh_.vertex(c))];
    }
    [[nodiscard]] int cornerPoint(int c) const { return mesh_.vertex(c); }
    [[nodiscard]] mesh::Vec3f position(int point) const { return points_[static_cast<std::size_t>(point)]; }
    [[nodiscard]] float doubleArea(int t) const { return doubleAreas_[static_cast<std::size_t>(t)]; }
    [[nodiscard]] int residualJump(int e) const { return residualJump_[static_cast<std::size_t>(e)]; }
    [[nodiscard]] bool cutEdge(int e) const { return cutEdges_.testBounded(static_cast<std::size_t>(e)); }
    // Whether the borders are frozen and `point` is a kept border point.
    [[nodiscard]] bool freezing() const noexcept { return frozen_ != nullptr; }
    [[nodiscard]] bool keptBorderPoint(int point) const
    {
        return frozen_ != nullptr && (*frozen_)[static_cast<std::size_t>(point)].value == 1;
    }
    // An edge with exactly one corner (false for -1).
    [[nodiscard]] bool borderEdge(int e) const { return e >= 0 && mesh_.isBoundaryEdge(e); }
    // A point with an edge that has exactly one corner.
    [[nodiscard]] bool borderPoint(int point) const;
    // Edge of the half edge at corner c (-1 for corners below 0 or without edge).
    [[nodiscard]] int cornerEdge(int c) const { return c >= 0 ? mesh_.cornerEdge(c) : -1; }
    // Whether the edges of the cut copy above edge e of the mesh are one (no seam).
    [[nodiscard]] bool seamless(int e) const;

    // Whether triangle t is turned over in (u, v): its signed area in (u, v)
    // below 0. Worked out once and kept in the trace.
    [[nodiscard]] bool flipped(int t);

    // Direction d given in the frame of triangle t2, in the frame of
    // triangle t1 (the triangle of `node`): unchanged across an edge of the
    // cut copy, turned by the quarter turns left across a shared edge of the
    // mesh, or summed over the edges around a shared point; then turned by
    // half a turn when exactly one of both is turned over. None when the
    // triangles share nothing or the point is not reached around.
    [[nodiscard]] std::optional<int> directionSlot(int t1, const TraceNode& node, int d, int t2);

    // Corner of `point` in triangle t (the first in corner order), or -1.
    [[nodiscard]] int cornerAt(int t, int point) const;

private:
    // Turns `s` by the quarter turns left across the edges of the fan around
    // the point of corner `from`, step by step, until corner `to` is
    // reached (true) or the fan ends (false, `s` as far as it got).
    [[nodiscard]] bool aroundPoint(int from, int to, int& s, bool forward) const;

    const mesh::Adjacency& mesh_;
    const field::CutCopy& copy_;
    const mesh::Adjacency& uvMesh_;
    std::span<const std::array<float, 2>> uv_;
    const std::vector<mesh::Vec3f>& points_;
    const std::vector<float>& doubleAreas_;
    std::span<const int> residualJump_;
    const mesh::BitArray& cutEdges_;
    const std::vector<mesh::FreezeMark>* frozen_;
    int triangles_;
    std::vector<std::uint8_t>& flipState_;
};

} // namespace areko::remesher::extract
