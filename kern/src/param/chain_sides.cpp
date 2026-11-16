#include "chain_sides.h"

#include <algorithm>
#include <array>

namespace areko::remesher::param {

using namespace areko::remesher::mesh;

namespace {

// From this many steps on a walk around a chain point must not return to its
// first corner, and from the second count on to no corner it has seen since.
constexpr int stepsBeforeReturnCheck = 300;
constexpr int stepsBeforeRepeatCheck = 500;

std::array<int, 2> cornerPair(const Adjacency& adjacency, int edge)
{
    const auto corners = adjacency.edgeCorners(edge);
    if (corners.size() == 1) {
        return {corners[0], -1};
    }
    if (corners.size() == 2) {
        return {corners[0], corners[1]};
    }
    return {-1, -1};
}

} // namespace

std::optional<ChainSides::CornerPairs> ChainSides::corners(const Adjacency& adjacency, std::span<const int> edges,
                                                           std::span<const int> points)
{
    CornerPairs out;
    if (edges.empty()) {
        return out;
    }
    const std::array<int, 2> first = cornerPair(adjacency, edges[0]);
    out.a.push_back(first[0]);
    out.b.push_back(first[1]);
    std::vector<int> seen;
    for (std::size_t k = 1; k < edges.size(); ++k) {
        const int pivot = points[k];
        const std::array<int, 2> pair = cornerPair(adjacency, edges[k]);
        const int begin = out.a.back();
        int corner = begin;
        int steps = 0;
        for (bool reached = false; !reached;) {
            if (corner < 0) {
                return std::nullopt;
            }
            // The other edge of the triangle at the chain point, and across it.
            const int across = adjacency.vertex(corner) == pivot ? adjacency.prev(corner) : adjacency.next(corner);
            const int edge = adjacency.cornerEdge(across);
            if (edge < 0) {
                return std::nullopt;
            }
            corner = adjacency.otherCorner(edge, across);
            if (corner < 0) {
                return std::nullopt;
            }
            ++steps;
            if (steps >= stepsBeforeReturnCheck) {
                if (corner == begin) {
                    return std::nullopt;
                }
                if (steps >= stepsBeforeRepeatCheck) {
                    if (steps == stepsBeforeRepeatCheck) {
                        seen.clear();
                    }
                    if (std::find(seen.begin(), seen.end(), corner) != seen.end()) {
                        return std::nullopt;
                    }
                    seen.push_back(corner);
                }
            }
            if (across == pair[1]) {
                out.a.push_back(pair[1]);
                out.b.push_back(pair[0]);
                reached = true;
            } else if (across == pair[0]) {
                out.a.push_back(pair[0]);
                out.b.push_back(pair[1]);
                reached = true;
            }
        }
    }
    return out;
}

ChainSides::Sides ChainSides::of(const Adjacency& adjacency, const field::CutCopy& copy,
                                 std::span<const int> residualJump, std::span<const int> edges,
                                 std::span<const int> points)
{
    Sides out;
    const std::optional<CornerPairs> pairs = corners(adjacency, edges, points);
    if (!pairs || edges.empty()) {
        return out;
    }
    out.found = true;
    out.jump = residualJump[static_cast<std::size_t>(edges[0])];
    if (pairs->a[0] != adjacency.edgeCorners(edges[0])[0]) {
        out.jump = -out.jump;
    }
    const Topology& copyCorners = copy.topology();
    const std::size_t count = points.size();
    for (std::size_t k = 0; k < count; ++k) {
        const std::size_t edge = k + 1 < count ? k : (count == 1 ? 0 : count - 2);
        const int a = pairs->a[edge];
        const int b = pairs->b[edge];
        if (a >= 0 && b >= 0) {
            auto copyPointAt = [&](int corner) {
                return copyCorners.vertex(adjacency.vertex(corner) == points[k] ? corner : adjacency.next(corner));
            };
            out.a.push_back(copyPointAt(a));
            out.b.push_back(copyPointAt(b));
        }
    }
    return out;
}

} // namespace areko::remesher::param
