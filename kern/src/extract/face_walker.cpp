#include "face_walker.h"

#include "seed_sectors.h"

#include <algorithm>

namespace areko::remesher::extract {

namespace {

// A walk is given up when it grows past this many nodes that moved.
constexpr int movedNodeLimit = 10;
// ... or past this many nodes in all.
constexpr int nodeLimit = 200;

} // namespace

bool FaceWalker::walk(std::array<float, 2> uv, int start, int direction, int triangle)
{
    TraceNodes& nodes = trace_.nodes;
    ++trace_.counts.walks;
    FaceLoop loop;
    loop.nodes.push_back(start);
    loop.directions.push_back(direction);
    loop.triangles.push_back(triangle);
    int unchanged = 0;
    const bool fromSeed = nodes[start].seed();
    bool firstStep = fromSeed;
    int d = direction;
    int t = triangle;
    for (;;) {
        const int k = loop.nodes.back();
        if (!nodes[k].seed()) {
            const auto slot = frame_.directionSlot(nodes[k].triangle, nodes[k], d, t);
            if (slot && nodes[k].done(*slot)) {
                return false;
            }
        }
        segments_.setCollecting(frame_.freezing());
        segments_.collected().clear();
        const std::array<float, 2> before = uv;
        const SegmentTracer::Result segment = segments_.run(k, d, t, uv, firstStep);
        if (segment.code != 0) {
            if (segment.code == 6) {
                trace_.deferredNodes.push_back(start);
                trace_.deferredDirections.push_back(direction);
            }
            return false;
        }
        const int k2 = segment.node;
        if (k2 < 0) {
            return false;
        }
        if (uv[0] == before[0] && uv[1] == before[1]) {
            ++unchanged;
        }
        if (k2 == k || (fromSeed && nodes[k2].seed() && nodes[k2].id < nodes[start].id)) {
            return false;
        }
        const auto m = static_cast<std::uint8_t>(segment.degenerate ? 1 : 0);
        if (segments_.collecting()) {
            for (const int f : segments_.collected()) {
                if (!(nodes[f].position == nodes[k].position) && !(nodes[f].position == nodes[k2].position)) {
                    loop.nodes.push_back(f);
                    loop.directions.push_back(d);
                    loop.triangles.push_back(nodes[f].triangle);
                    loop.degenerate.push_back(m);
                }
            }
        }
        if (loop.nodes.front() == k2) {
            if (loop.nodes.size() < 3) {
                return false;
            }
            if (fromSeed && nodes[loop.nodes[1]].seed() && faceExists(trace_, loop)) {
                return false;
            }
            loop.degenerate.push_back(m);
            settle(loop);
            closer_.close(loop, trace_);
            ++trace_.counts.closedWalks;
            return true;
        }
        const TraceNode& reached = nodes[k2];
        if (!(reached.seed() && reached.kind == node_kind::point)) {
            d = (segment.arrival + (frame_.flipped(t) ? 3 : 1)) % 4;
        } else {
            d = segment.arrival;
            if (!SeedSectors::choose(frame_, reached.reference, t, d, uv) || t == -1) {
                return false;
            }
            const auto corners = frame_.mesh().pointCorners(reached.reference);
            if (std::none_of(corners.begin(), corners.end(),
                             [&](int c) { return frame_.mesh().cornerFace(c) == t; })) {
                return false;
            }
        }
        loop.nodes.push_back(k2);
        loop.directions.push_back(d);
        loop.triangles.push_back(t);
        loop.degenerate.push_back(m);
        if (tooLong(loop, unchanged)) {
            ++trace_.counts.tooLong;
            return false;
        }
        firstStep = false;
    }
}

bool FaceWalker::tooLong(const FaceLoop& loop, int unchanged)
{
    const int n = static_cast<int>(loop.nodes.size());
    if (n <= movedNodeLimit) {
        return false;
    }
    if (n - unchanged > movedNodeLimit) {
        return true;
    }
    // A node met a second time: looked for among the nodes before the one
    // just added.
    const auto earlier = loop.nodes.end() - 1;
    if (std::find(loop.nodes.begin(), earlier, loop.nodes.back()) != earlier) {
        return true;
    }
    return n > nodeLimit;
}

bool FaceWalker::faceExists(const Trace& trace, const FaceLoop& loop)
{
    const mesh::Csr& polygons = trace.polygons;
    const std::size_t n = loop.nodes.size();
    for (int r = 0; r < polygons.rows(); ++r) {
        const auto polygon = polygons.row(r);
        if (polygon.size() != n || n == 0) {
            continue;
        }
        bool same = true;
        for (std::size_t i = 0; i < n && same; ++i) {
            const TraceNode& node = trace.nodes[loop.nodes[i]];
            if (!node.seed()) {
                return false;
            }
            same = node.id == polygon[i];
        }
        if (same) {
            return true;
        }
    }
    return false;
}

void FaceWalker::settle(const FaceLoop& loop)
{
    TraceNodes& nodes = trace_.nodes;
    for (std::size_t i = 0; i < loop.nodes.size(); ++i) {
        const int k = loop.nodes[i];
        if (nodes[k].id == -1) {
            nodes[k].id = nodes.takeId();
        }
        if (!nodes[k].seed()) {
            if (const auto slot = frame_.directionSlot(nodes[k].triangle, nodes[k], loop.directions[i], loop.triangles[i])) {
                nodes[k].flags = static_cast<std::uint8_t>(nodes[k].flags | (1u << *slot));
            }
        }
    }
}

} // namespace areko::remesher::extract
