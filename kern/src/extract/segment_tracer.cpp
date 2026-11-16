#include "segment_tracer.h"

#include "node_search.h"
#include "util/rounding.h"

#include <cfloat>
#include <cmath>
#include <utility>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

using Uv = std::array<float, 2>;

// Widened tolerance for the way out of a triangle: the float 5e-5 in double.
constexpr double outTolerance = static_cast<double>(5e-5f);

bool unusable(float v) { return std::isnan(v) || std::fabs(v) > FLT_MAX; }

bool sameUv(const Uv& p, const Uv& q) { return p[0] == q[0] && p[1] == q[1]; }

// Distance in (u, v) between two points, in float.
float uvDistance(const Uv& p, const Uv& q)
{
    const float du = p[0] - q[0];
    const float dv = p[1] - q[1];
    return std::sqrt(du * du + dv * dv);
}

// a + (b - a) * f per coordinate, in float.
Vec3f along(Vec3f a, Vec3f b, float f) { return {(b.x - a.x) * f + a.x, (b.y - a.y) * f + a.y, (b.z - a.z) * f + a.z}; }

std::size_t at(int k) { return static_cast<std::size_t>(k); }

} // namespace

struct SegmentTracer::Walk {
    Walk(const TraceNode& from, int s, int d, int& triangle, Uv& coordinates, bool first)
        : start(s), direction(d), tri(triangle), uv(coordinates), firstStep(first), startSeed(from.seed()),
          startKind(from.kind), startPoint(from.reference), uvStart(coordinates), triStart(triangle), a(d % 2),
          c(1 - a), arrival(d), curF(coordinates[at(a)]), cur(curF), basis(cur), upper(cur + 1.0), lower(cur - 1.0)
    {
    }

    int start;
    int direction;
    int& tri;
    Uv& uv;
    bool firstStep;
    bool startSeed;
    int startKind;
    int startPoint;
    Uv uvStart;
    int triStart;
    int counter = 0;
    // Parameters of the last way out found upwards and downwards; they are
    // kept over the whole segment.
    double maxT = 0.0;
    double minT = 0.0;
    bool overCorner = false;
    int a; // moving coordinate
    int c; // fixed coordinate
    int arrival;
    int entry = -1;      // corner the segment entered the triangle over
    int prevTri = -1;    // triangle left last
    int hPrev = -1;      // half edge or corner the triangle was entered over
    int cornerMark = -1;
    int h = -1;          // corner the line leaves the triangle through
    int side = -1;       // edge of the way out (0: 01, 1: 12, 2: 20)
    float curF;
    double cur;
    double basis;
    double upper;
    double lower;
};

struct SegmentTracer::Cell {
    int base = 0;
    std::array<Uv, 3> p{};
    float fixed = 0.0f;
    // Parameter of the crossing of the line with edge k (from corner k).
    std::array<double, 3> t{};
    double x = 0.0;    // moving coordinate at the way out
    double tOut = 0.0; // edge parameter of the way out
    int corner = -1;   // corner the line leaves through
    bool followsEdge = false;

    [[nodiscard]] float coordinate(int k, int axis) const { return p[at(k)][at(axis)]; }
};

struct SegmentTracer::Target {
    float value = 0.0f;
    int result = 0; // 0 the start value again, 1 one above, -1 one below
    bool exact = false;
};

SegmentTracer::Result SegmentTracer::run(int start, int direction, int& triangle, Uv& uv, bool firstStep)
{
    // A segment already followed from this node in this direction.
    const TraceNode& node = trace_.nodes[start];
    if (!node.seed()) {
        if (const auto slot = frame_.directionSlot(node.triangle, node, direction, triangle)) {
            const NodeLink& link = node.links[static_cast<std::size_t>(*slot)];
            if (link.node >= 0) {
                triangle = link.triangle;
                uv = link.uv;
                ++trace_.counts.linkedSegments;
                return {0, link.node, link.direction, link.degenerate};
            }
        }
    }
    const Result result = follow(start, direction, triangle, uv, firstStep);
    ++trace_.counts.segments[static_cast<std::size_t>(result.code)];
    return result;
}

SegmentTracer::Result SegmentTracer::follow(int s, int d, int& tri, Uv& uv, bool firstStep)
{
    Walk walk(trace_.nodes[s], s, d, tri, uv, firstStep);
    for (;;) {
        ++trace_.segmentVisits[static_cast<std::size_t>(walk.tri)];
        Cell cell = enter(walk);
        if (const auto stop = wayOut(walk, cell)) {
            return *stop;
        }
        if (const auto stop = throughCorner(walk, cell)) {
            return *stop;
        }
        if (const auto target = reached(walk, cell)) {
            return finish(walk, cell, *target);
        }
        if (cell.x < walk.lower || walk.upper < cell.x) {
            return {2, -1, walk.arrival, false};
        }
        collectBorderPoint(walk);
        ++walk.counter;
        if (const auto stop = crossOver(walk, cell)) {
            return *stop;
        }
        walk.curF = walk.uv[at(walk.a)];
        walk.cur = walk.curF;
        walk.hPrev = walk.h;
    }
}

SegmentTracer::Cell SegmentTracer::enter(const Walk& walk) const
{
    Cell cell;
    cell.base = 3 * walk.tri;
    cell.p = {frame_.cornerUv(cell.base), frame_.cornerUv(cell.base + 1), frame_.cornerUv(cell.base + 2)};
    cell.fixed = walk.uv[at(walk.c)];
    for (int k = 0; k < 3; ++k) {
        const int n = (k + 1) % 3;
        cell.t[at(k)] = (wide(cell.fixed) - wide(cell.coordinate(k, walk.c))) /
                        (wide(cell.coordinate(n, walk.c)) - wide(cell.coordinate(k, walk.c)));
    }
    return cell;
}

std::optional<SegmentTracer::Result> SegmentTracer::wayOut(Walk& walk, Cell& cell) const
{
    const int a = walk.a;
    double tolMax = 0.0;
    double tolMin = 0.0;
    int attempt = 0;
    bool atStart = false;
    double maxX = walk.cur;
    double minX = walk.cur;
    int maxJ = -1;
    int minJ = -1;
    for (;;) {
        atStart = sameUv(walk.uv, walk.uvStart) && !walk.startSeed;
        std::array<bool, 3> blocked{};
        if (!atStart) {
            if (walk.entry >= 0) {
                blocked[at(walk.entry % 3)] = true;
            } else if (walk.startSeed && walk.startKind == node_kind::point) {
                const int i = frame_.cornerAt(walk.tri, walk.startPoint);
                if (i < 0) {
                    return Result{0, -1, walk.arrival, false};
                }
                blocked[at(i % 3)] = true;
                blocked[at((i + 2) % 3)] = true;
            } else if (walk.overCorner && walk.hPrev >= 0) {
                // Never taken: in the first triangle the line has not left
                // over a corner yet before its way out is chosen.
                blocked[at(walk.hPrev % 3)] = true;
                blocked[at((walk.hPrev + 2) % 3)] = true;
            }
        }
        maxX = walk.cur;
        minX = walk.cur;
        maxJ = -1;
        minJ = -1;
        for (int k = 0; k < 3; ++k) {
            const double tk = cell.t[at(k)];
            if (!(tk >= 0.0 && tk <= 1.0) || blocked[at(k)]) {
                continue;
            }
            const int n = (k + 1) % 3;
            const double x = wide(cell.coordinate(k, a)) + wide(cell.coordinate(n, a) - cell.coordinate(k, a)) * tk;
            if (x >= maxX - tolMax) {
                maxX = x;
                walk.maxT = tk;
                maxJ = k;
            }
            if (x <= minX + tolMin) {
                minX = x;
                walk.minT = tk;
                minJ = k;
            }
        }
        if (maxJ >= 0 || minJ >= 0) {
            break;
        }
        if (attempt == 0) {
            (walk.arrival < 2 ? tolMin : tolMax) = outTolerance;
            attempt = 1;
        } else if (attempt == 1) {
            tolMax = outTolerance;
            tolMin = outTolerance;
            attempt = 2;
        } else if (attempt == 2) {
            // Parameters anew from the moving coordinate where it is close;
            // a gap that is not a number counts as far.
            std::array<double, 3> fresh{};
            std::array<double, 3> gap{};
            for (int k = 0; k < 3; ++k) {
                const int n = (k + 1) % 3;
                fresh[at(k)] = wide((walk.curF - cell.coordinate(k, a)) / (cell.coordinate(n, a) - cell.coordinate(k, a)));
                gap[at(k)] = std::fabs(cell.t[at(k)] - fresh[at(k)]);
            }
            if (!(gap[0] < 0.01) && !(gap[1] < 0.01) && !(gap[2] < 0.01)) {
                return Result{4, -1, walk.arrival, false};
            }
            for (std::size_t k = 0; k < 3; ++k) {
                if (gap[k] < 0.01) {
                    cell.t[k] = fresh[k];
                }
            }
            attempt = 3;
        } else {
            return Result{4, -1, walk.arrival, false};
        }
    }

    const bool upwards = atStart ? walk.arrival < 2
                                 : !(std::fabs(static_cast<float>(walk.cur - minX)) >
                                     std::fabs(static_cast<float>(walk.cur - maxX))) ||
                                       minJ == -1;
    if (upwards) {
        cell.tOut = walk.maxT;
        cell.x = maxX;
        walk.side = maxJ;
    } else {
        cell.tOut = walk.minT;
        cell.x = minX;
        walk.side = minJ;
    }
    if (atStart && walk.side == -1) {
        return Result{4, -1, walk.arrival, false};
    }
    return std::nullopt;
}

std::optional<SegmentTracer::Result> SegmentTracer::throughCorner(Walk& walk, Cell& cell) const
{
    const int a = walk.a;
    const int c = walk.c;
    int corner = -1;
    if (cell.tOut == 0.0) {
        corner = walk.side;
    } else if (cell.tOut == 1.0) {
        corner = (walk.side + 1) % 3;
    } else {
        for (int k = 0; k < 3 && corner < 0; ++k) {
            corner = wide(cell.coordinate(k, a)) == cell.x && cell.coordinate(k, c) == cell.fixed ? k : -1;
        }
    }
    cell.corner = corner;
    if (corner >= 0) {
        if (walk.firstStep && walk.startKind == node_kind::point &&
            frame_.cornerPoint(cell.base + corner) == walk.startPoint) {
            return Result{5, -1, walk.arrival, false};
        }
        walk.h = cell.base + corner;
        if (walk.overCorner) {
            if (walk.cornerMark / 3 == walk.tri) {
                return Result{3, -1, walk.arrival, false};
            }
            if (walk.h != walk.hPrev) {
                walk.cornerMark = walk.h;
            }
        } else {
            cell.followsEdge = cell.fixed == frame_.cornerUv(frame_.mesh().next(walk.h))[at(c)];
            walk.overCorner = true;
            walk.cornerMark = walk.h;
        }
    } else {
        if (walk.firstStep && cell.x == walk.cur) {
            return Result{5, -1, walk.arrival, false};
        }
        walk.cornerMark = -1;
        walk.overCorner = false;
        walk.h = walk.hPrev;
    }
    return std::nullopt;
}

std::optional<SegmentTracer::Target> SegmentTracer::reached(const Walk& walk, const Cell& cell)
{
    const double x = cell.x;
    auto past = [&](float z) { return (walk.curF < z && wide(z) <= x) || (z < walk.curF && x <= wide(z)); };
    const float onBase = static_cast<float>(walk.basis);
    const float onUpper = static_cast<float>(walk.upper);
    const float onLower = static_cast<float>(walk.lower);
    const auto& p = cell.p;
    if (past(onBase)) {
        return Target{onBase, 0, false};
    }
    if (onBase == walk.curF && wide(onBase) == x && (sameUv(p[0], p[1]) || sameUv(p[1], p[2]) || sameUv(p[0], p[2]))) {
        return Target{onBase, 0, true};
    }
    if (past(onUpper)) {
        return Target{onUpper, 1, false};
    }
    if (x == wide(onUpper)) {
        return Target{onUpper, 1, true};
    }
    if (past(onLower)) {
        return Target{onLower, -1, false};
    }
    if (x == wide(onLower)) {
        return Target{onLower, -1, true};
    }
    return std::nullopt;
}

SegmentTracer::Result SegmentTracer::finish(Walk& walk, const Cell& cell, const Target& target)
{
    TraceNodes& nodes = trace_.nodes;
    const double x = cell.x;
    const bool extended = target.exact || x == wide(target.value) || cell.corner >= 0;
    walk.uv[at(walk.a)] = target.value;
    Vec3f position;
    if (!place(frame_, walk.tri, cell.corner, walk.uv, position)) {
        return {5, -1, walk.arrival, false};
    }
    // A value that is not a number counts as not above.
    if (!(x > walk.cur)) {
        if (x != walk.cur && walk.arrival < 2) {
            walk.arrival += 2;
        }
    } else if (walk.arrival > 1) {
        walk.arrival -= 2;
    }
    const bool degenerate = target.result == 0;
    if (degenerate && stackPhase_) {
        return {6, -1, walk.arrival, true};
    }
    int end = NodeSearch::find(frame_, trace_, trace_.nodeToleranceSq, walk.tri, walk.side, position, extended);
    if (end >= 0) {
        const TraceNode& found = nodes[end];
        if (!found.seed() && found.triangle != walk.tri) {
            if (const auto slot = frame_.directionSlot(found.triangle, found, walk.arrival, walk.tri)) {
                walk.tri = found.triangle;
                walk.uv = found.uv;
                walk.arrival = *slot;
            }
        }
    } else {
        end = newNode(walk, cell, position);
    }
    // Not while border points of a frozen border were collected on it.
    if (!(collecting_ && !collected_.empty())) {
        connect(walk, end, degenerate);
    }
    return {0, end, walk.arrival, degenerate};
}

int SegmentTracer::newNode(const Walk& walk, const Cell& cell, Vec3f position)
{
    TraceNodes& nodes = trace_.nodes;
    const int k = nodes.create(walk.tri, position, walk.uv, true);
    TraceNode& made = nodes[k];
    if (walk.overCorner) {
        if (sameUv(walk.uv, frame_.cornerUv(walk.h))) {
            made.kind = node_kind::point;
            made.reference = frame_.cornerPoint(walk.h);
            if (frame_.borderPoint(made.reference)) {
                made.flags |= node_flag::border;
            }
        } else if (walk.hPrev != -1) {
            const int pa = frame_.cornerPoint(walk.h);
            const int pb = frame_.cornerPoint(walk.hPrev);
            if (pa != pb && frame_.borderEdge(frame_.mesh().edgeBetween(pa, pb))) {
                made.flags |= node_flag::border;
            }
        }
    }
    if (cell.x == wide(walk.uv[at(walk.a)]) && frame_.borderEdge(frame_.cornerEdge(cell.base + walk.side))) {
        made.flags |= node_flag::border;
    }
    NodeSearch::snapToEdge(frame_, made);
    if ((nodes[walk.start].flags & node_flag::border) != 0 && (made.flags & node_flag::border) == 0) {
        borderAlongLine(k, walk.a);
    }
    return k;
}

void SegmentTracer::connect(const Walk& walk, int end, bool degenerate)
{
    TraceNodes& nodes = trace_.nodes;
    if (!walk.startSeed) {
        const TraceNode& from = nodes[walk.start];
        if (const auto slot = frame_.directionSlot(from.triangle, from, walk.direction, walk.triStart)) {
            NodeLink& link = nodes[walk.start].links[static_cast<std::size_t>(*slot)];
            if (link.node < 0) {
                link = {end, walk.arrival, walk.tri, walk.uv, degenerate};
            }
        }
    }
    const TraceNode& to = nodes[end];
    if (!to.seed()) {
        if (const auto slot = frame_.directionSlot(to.triangle, to, (walk.arrival + 2) % 4, walk.tri)) {
            NodeLink& link = nodes[end].links[static_cast<std::size_t>(*slot)];
            if (link.node < 0) {
                link = {walk.start, (walk.direction + 2) % 4, walk.triStart, walk.uvStart, degenerate};
            }
        }
    }
}

void SegmentTracer::collectBorderPoint(const Walk& walk)
{
    if (!collecting_ || !walk.overCorner) {
        return;
    }
    TraceNodes& nodes = trace_.nodes;
    const int q = frame_.cornerPoint(walk.h);
    const bool sameAsLast = !collected_.empty() && nodes[collected_.back()].kind == node_kind::point &&
                            nodes[collected_.back()].reference == q;
    if (frame_.keptBorderPoint(q) && !sameAsLast && !(walk.startKind == node_kind::point && q == walk.startPoint)) {
        const int f = nodes.create(walk.tri, frame_.position(q), frame_.cornerUv(walk.h), false);
        nodes[f].kind = node_kind::point;
        nodes[f].reference = q;
        nodes[f].flags = node_flag::border | node_flag::frozenBorder;
        collected_.push_back(f);
    }
}

std::optional<SegmentTracer::Result> SegmentTracer::crossOver(Walk& walk, Cell& cell) const
{
    const Adjacency& mesh = frame_.mesh();
    do {
        int edge = -1;
        if (walk.overCorner) {
            edge = frame_.cornerEdge(walk.h);
            const int across = edge >= 0 ? mesh.otherCorner(edge, walk.h) : -1;
            cell.tOut = 0.0;
            walk.hPrev = walk.h;
            walk.h = across >= 0 ? mesh.next(across) : -1;
        } else {
            walk.hPrev = cell.base + walk.side;
            edge = frame_.cornerEdge(walk.hPrev);
        }
        const int other = edge >= 0 ? mesh.otherCorner(edge, walk.hPrev) : -1;
        if (other < 0) {
            trace_.nodes[walk.start].flags |= node_flag::ranIntoBorder;
            return Result{1, -1, walk.arrival, false};
        }
        if (other / 3 == walk.prevTri) {
            return Result{3, -1, walk.arrival, false};
        }
        walk.prevTri = walk.tri;
        walk.tri = other / 3;
        walk.entry = other;
        if (frame_.seamless(edge)) {
            walk.uv[at(walk.a)] = static_cast<float>(cell.x);
            if (walk.counter == 5001) {
                return Result{3, -1, walk.arrival, false};
            }
        } else {
            if (walk.counter > 1000) {
                return Result{3, -1, walk.arrival, false};
            }
            const Uv& q1 = frame_.cornerUv(other);
            const Uv& q2 = frame_.cornerUv(mesh.next(other));
            walk.uv[0] = static_cast<float>(wide(q2[0]) * (1.0 - cell.tOut) + wide(q1[0]) * cell.tOut);
            walk.uv[1] = static_cast<float>(wide(q2[1]) * (1.0 - cell.tOut) + wide(q1[1]) * cell.tOut);
            int r = frame_.residualJump(edge);
            if (walk.prevTri != mesh.edgeCorners(edge)[0] / 3) {
                r = -r;
            }
            const int before = walk.arrival;
            walk.arrival = mod4(before - r);
            walk.a = walk.arrival % 2;
            walk.c = 1 - walk.a;
            const float shift = static_cast<float>(cell.x - walk.basis);
            const float now = walk.uv[at(walk.a)];
            const float moved = (before < 2) != (walk.arrival < 2) ? shift + now : now - shift;
            cell.x = wide(now);
            // The new start value is rounded in double, the fixed coordinate in float.
            walk.basis = static_cast<double>(util::roundHalfAway(static_cast<double>(moved)));
            walk.upper = walk.basis + 1.0;
            walk.lower = walk.basis - 1.0;
            walk.uv[at(walk.c)] = static_cast<float>(util::roundHalfAway(walk.uv[at(walk.c)]));
        }
    } while (std::exchange(cell.followsEdge, false));
    return std::nullopt;
}

bool SegmentTracer::place(const TraceFrame& frame, int t, int corner, const Uv& uv, Vec3f& out)
{
    const int base = 3 * t;
    const Uv& p0 = frame.cornerUv(base);
    const Uv& p1 = frame.cornerUv(base + 1);
    const Uv& p2 = frame.cornerUv(base + 2);
    const Vec3f q0 = frame.position(frame.cornerPoint(base));
    const Vec3f q1 = frame.position(frame.cornerPoint(base + 1));
    const Vec3f q2 = frame.position(frame.cornerPoint(base + 2));

    // Barycentric in (u, v): differences in float, the rest in double.
    const double du1 = wide(p1[0] - p0[0]);
    const double dv1 = wide(p1[1] - p0[1]);
    const double du2 = wide(p2[0] - p0[0]);
    const double dv2 = wide(p2[1] - p0[1]);
    const double du = wide(uv[0] - p0[0]);
    const double dv = wide(uv[1] - p0[1]);
    const double det = du1 * dv2 - dv1 * du2;
    const double l1 = (dv2 / det) * du + (-du2 / det) * dv;
    const double l2 = (-dv1 / det) * du + (du1 / det) * dv;
    auto mix = [&](float a0, float a1, float a2) {
        return static_cast<float>((wide(a1 - a0) * l1 + wide(a0)) + wide(a2 - a0) * l2);
    };
    out = {mix(q0.x, q1.x, q2.x), mix(q0.y, q1.y, q2.y), mix(q0.z, q1.z, q2.z)};
    if (!unusable(out.x)) {
        return true;
    }
    if (corner >= 0) {
        out = frame.position(frame.cornerPoint(base + corner));
        return true;
    }
    if (sameUv(p1, p0) && sameUv(p1, p2)) {
        out = {(q0.x + q1.x + q2.x) / 3.0f, (q0.y + q1.y + q2.y) / 3.0f, (q0.z + q1.z + q2.z) / 3.0f};
    } else {
        // Along the longest side of the triangle in (u, v); a length that is
        // not a number counts as not longer.
        const float l01 = uvDistance(p0, p1);
        const float l12 = uvDistance(p1, p2);
        const float l20 = uvDistance(p2, p0);
        if (!(l12 < l01)) {
            if (!(l20 < l12)) {
                out = along(q0, q2, uvDistance(uv, p0) / l20);
            } else {
                out = along(q1, q2, uvDistance(uv, p1) / l12);
            }
        } else if (!(l20 < l01)) {
            out = along(q0, q2, uvDistance(uv, p0) / l20);
        } else {
            out = along(q0, q1, uvDistance(uv, p0) / l01);
        }
    }
    return !unusable(out.x) && !unusable(out.y) && !unusable(out.z);
}

void SegmentTracer::borderAlongLine(int node, int a)
{
    TraceNode& n = trace_.nodes[node];
    if (n.triangle < 0 || n.triangle >= frame_.triangleCount()) {
        return;
    }
    const int base = 3 * n.triangle;
    const int fixedAxis = 1 - a;
    for (int k = 0; k < 3; ++k) {
        const Uv& pk = frame_.cornerUv(base + k);
        const Uv& pn = frame_.cornerUv(base + (k + 1) % 3);
        const float value = n.uv[static_cast<std::size_t>(fixedAxis)];
        const bool onLine = value == pk[static_cast<std::size_t>(fixedAxis)] &&
                            value == pn[static_cast<std::size_t>(fixedAxis)];
        if (onLine && frame_.borderEdge(frame_.cornerEdge(base + k))) {
            n.flags |= node_flag::border;
        }
    }
}

} // namespace areko::remesher::extract
