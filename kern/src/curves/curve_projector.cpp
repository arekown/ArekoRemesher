#include "curve_projector.h"

#include "prep/polygon_triangulator.h"
#include "util/introsort.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <deque>
#include <numeric>

namespace areko::remesher::curves {

using namespace areko::remesher::mesh;

namespace {

constexpr int maxWalkSteps = 1001;
constexpr float absoluteTolerance = 1e-5f;


// A located point: its position and where it lies.
struct Spot {
    Vec3f position;
    Location location;
};

// Result of the closest approach of a line (c, direction) to an edge line.
struct Approach {
    double s = 0.0; // along the direction
    double t = 0.0; // along the edge, 0 first and 1 second end point
    bool parallel = false;
    double distance = 0.0; // between the lines when parallel
};

// One crossing candidate of the straight walk: position s along the
// direction, parameter t on the edge, the face crossed and the point
// (1 - t) p0 + t p1 in double.
struct Crossing {
    bool valid = false;
    double s = -1e6;
    double t = 0.0;
    int edge = -1;
    int face = -1;
    Vec3d point;
};

// Unit normal of every face of a polygon mesh: Newell sum in double, its
// squared length rounded to float and its root taken in float; above 1e-26
// every coordinate is divided in double and rounded to float, otherwise the
// sum stays as it is (rounded to float).
std::vector<Vec3f> polygonNormals(const Mesh& mesh)
{
    const Topology& topology = mesh.topology();
    std::vector<Vec3f> out(static_cast<std::size_t>(topology.faceCount()));
    std::vector<Vec3f> corners;
    for (int f = 0; f < topology.faceCount(); ++f) {
        corners.clear();
        for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
            corners.push_back(mesh.points()[static_cast<std::size_t>(topology.vertex(c))]);
        }
        Vec3d n = prep::PolygonTriangulator::newellSum(corners);
        const auto squared = static_cast<float>((n.x * n.x + n.y * n.y) + n.z * n.z);
        const auto len = static_cast<double>(std::sqrt(squared));
        if (len > 1e-26) {
            n = {wide(static_cast<float>(n.x / len)), wide(static_cast<float>(n.y / len)),
                 wide(static_cast<float>(n.z / len))};
        }
        out[static_cast<std::size_t>(f)] = {static_cast<float>(n.x), static_cast<float>(n.y), static_cast<float>(n.z)};
    }
    return out;
}

class Tracer {
public:
    // Without `walks` the tracer only looks for paths over edges.
    Tracer(const Mesh& mesh, const CurveProjector::Options& options, bool walks = true)
        : mesh_(mesh), adjacency_(mesh.adjacency()), topology_(mesh.topology()), points_(mesh.points()),
          options_(options)
    {
        if (walks) {
            if (topology_.mode() == FaceMode::Triangles) {
                normals_ = &mesh.triangleGeometry().normals();
            } else {
                ownNormals_ = polygonNormals(mesh);
                normals_ = &ownNormals_;
            }
            tolerance_ = options.tolerance > 0.0f ? options.tolerance
                                                  : std::max(diagonal(mesh.bounds()) * 1e-6f, 1e-8f);
            tolerance2_ = tolerance_ * tolerance_;
        }
    }

    [[nodiscard]] std::vector<CurveSegment> project(const LocatedPolyline& curve) const;
    [[nodiscard]] bool edgePath(const Spot& start, const Spot& target, CurveSegment& out) const;

private:
    [[nodiscard]] Vec3f at(int v) const { return points_[static_cast<std::size_t>(v)]; }
    // Corner `steps` places after corner c in its face (cyclic).
    [[nodiscard]] int cornerAfter(int c, int steps) const
    {
        const int f = adjacency_.cornerFace(c);
        const int begin = topology_.faceBegin(f);
        return begin + (c - begin + steps) % topology_.faceSize(f);
    }
    [[nodiscard]] int firstEdge(int v) const
    {
        const auto edges = adjacency_.pointEdges(v);
        return edges.empty() ? -1 : edges.front();
    }
    // Parameter of end point v on edge e: 0 for the first end point, else 1.
    [[nodiscard]] float endParameter(int e, int v) const { return adjacency_.edgeVertices(e)[0] == v ? 0.0f : 1.0f; }

    void locate(Spot& spot) const;
    [[nodiscard]] bool connect(const Spot& from, const Spot& to, CurveSegment& out) const;
    [[nodiscard]] bool walk(const Spot& start, const Spot& target, CurveSegment& out, bool& backward) const;
    void cornerFallback(const Location& here, const Vec3d& c, const Vec3d& direction, const std::vector<int>& halfEdges,
                        Crossing& best) const;
    [[nodiscard]] bool candidates(const Location& here, const Spot& target, std::vector<int>& halfEdges,
                                  CurveSegment& out) const;
    [[nodiscard]] Approach approach(const Vec3d& c, const Vec3d& direction, Vec3f a, Vec3f b) const;
    [[nodiscard]] float pathLength(const Spot& start, const CurveSegment& path, const Spot& target) const;
    [[nodiscard]] int commonFace(int e0, int e1) const;
    [[nodiscard]] float chordLimit(const Spot& start, const Spot& target) const;
    [[nodiscard]] bool layOntoChord(const Spot& start, const Spot& target, CurveSegment& path, bool giveUp) const;

    // Entry of a path of the walk: every entry names the face crossed to reach it.
    static void appendPathEntry(CurveSegment& out, int face, int edge, float t, Vec3f p)
    {
        out.faces.push_back(face);
        out.edges.push_back(edge);
        out.points.push_back(p);
        out.t.push_back(t);
    }
    static void appendEntry(CurveSegment& out, int face, int edge, float t, Vec3f p)
    {
        if (!out.empty()) {
            out.faces.push_back(face);
        }
        out.edges.push_back(edge);
        out.points.push_back(p);
        out.t.push_back(t);
    }
    // Same as the last entry (edge and parameter): nothing to add.
    static bool repeatsLast(const CurveSegment& out, int edge, float t)
    {
        return !out.empty() && out.edges.back() == edge && static_cast<double>(out.t.back()) == static_cast<double>(t);
    }

    const Mesh& mesh_;
    const Adjacency& adjacency_;
    const Topology& topology_;
    const std::vector<Vec3f>& points_;
    CurveProjector::Options options_;
    std::vector<Vec3f> ownNormals_;
    const std::vector<Vec3f>* normals_ = nullptr;
    float tolerance_ = 0.0f;
    float tolerance2_ = 0.0f;
};

void Tracer::locate(Spot& spot) const
{
    if (spot.location.kind != Location::Face) {
        return;
    }
    const int f = spot.location.index;
    const Vec3f p = spot.position;
    for (int c = 3 * f; c < 3 * f + 3; ++c) {
        const int v = topology_.vertex(c);
        const Vec3f d = p - at(v);
        if (dot(d, d) < tolerance2_) {
            spot.location = {Location::Point, v};
        }
    }
    if (spot.location.kind == Location::Point) {
        spot.position = at(spot.location.index);
        return;
    }
    Vec3f foot;
    for (int c = 3 * f; c < 3 * f + 3; ++c) {
        const int e = adjacency_.cornerEdge(c);
        if (e < 0) {
            continue;
        }
        const auto ends = adjacency_.edgeVertices(e);
        const Vec3f a = at(ends[0]);
        Vec3f dir = at(ends[1]) - a;
        const float len = length(dir);
        if (len <= -1e-14f || 1e-14f <= len) {
            dir = {dir.x / len, dir.y / len, dir.z / len};
        } else {
            dir = {};
        }
        const Vec3f w = p - a;
        const float along = dot(w, dir);
        const Vec3f off = w - dir * along;
        if (dot(off, off) < tolerance2_) {
            spot.location = {Location::Edge, e};
            foot = p - off;
        }
    }
    if (spot.location.kind != Location::Edge) {
        return;
    }
    spot.position = foot;
    // Onto an end point of that edge when it lies within the tolerance.
    const auto ends = adjacency_.edgeVertices(spot.location.index);
    for (const int v : ends) {
        const Vec3f d = foot - at(v);
        if (dot(d, d) < tolerance2_) {
            spot.location = {Location::Point, v};
        }
    }
    if (spot.location.kind == Location::Point) {
        spot.position = at(spot.location.index);
    }
}

// Face shared by the edges e0 and e1. Two inner edges: the face of the first
// corner of e1 when one of the faces of e0 is it, else that of its second
// corner; otherwise the first match with the corners of e0 outside.
int Tracer::commonFace(int e0, int e1) const
{
    const auto a = adjacency_.edgeCorners(e0);
    const auto b = adjacency_.edgeCorners(e1);
    if (a.size() == 2 && b.size() == 2) {
        const int a0 = adjacency_.cornerFace(a[0]);
        const int a1 = adjacency_.cornerFace(a[1]);
        for (const int c : b) {
            const int f = adjacency_.cornerFace(c);
            if (a0 == f || a1 == f) {
                return f;
            }
        }
        return -1;
    }
    for (const int c : adjacency_.edgeCorners(e0)) {
        const int f = adjacency_.cornerFace(c);
        for (const int d : adjacency_.edgeCorners(e1)) {
            if (adjacency_.cornerFace(d) == f) {
                return f;
            }
        }
    }
    return -1;
}

bool Tracer::connect(const Spot& from, const Spot& to, CurveSegment& out) const
{
    const int kf = from.location.kind;
    const int kt = to.location.kind;
    const int a = from.location.index;
    const int b = to.location.index;
    auto add = [&](int face, int edge, float t, Vec3f p) {
        if (!repeatsLast(out, edge, t)) {
            appendEntry(out, face, edge, t, p);
        }
        return true;
    };
    if (kf == Location::Point && kt == Location::Point) {
        const int e = adjacency_.edgeBetween(b, a);
        if (e >= 0) {
            return add(-1, e, endParameter(e, b), at(b));
        }
        return false;
    }
    if (kf == Location::Face && kt == Location::Face) {
        if (a == b) {
            return true; // face points are not kept
        }
        const int e = adjacency_.commonEdge(b, a);
        if (e < 0) {
            return false;
        }
        // Where the line between the two points crosses the shared edge.
        const auto ends = adjacency_.edgeVertices(e);
        const Vec3f p0 = at(ends[0]);
        const Vec3f p1 = at(ends[1]);
        Vec3d dir = wide(to.position) - wide(from.position);
        normalizeInPlace(dir);
        const Approach hit = approach(wide(from.position), dir, p0, p1);
        float t = 0.0f;
        if (!hit.parallel) {
            t = static_cast<float>(hit.t);
        } else {
            if (hit.distance > 1e-4) {
                return false;
            }
            const Vec3d mid = wide((p0 + p1) * 0.5f) - wide(p0);
            const double edgeLength = length(wide(p1) - wide(p0));
            t = static_cast<float>(std::min(length(mid) / edgeLength, 1.0));
        }
        const auto tf = t;
        const auto tr = static_cast<float>(1.0 - static_cast<double>(t));
        const Vec3f p{p1.x * tf + p0.x * tr, p1.y * tf + p0.y * tr, p1.z * tf + p0.z * tr};
        return add(a, e, t, p);
    }
    if (kf == Location::Face && kt == Location::Point) {
        for (int c = 3 * a; c < 3 * a + 3; ++c) {
            if (topology_.vertex(c) == b) {
                const int e = adjacency_.cornerEdge(c);
                return add(a, e, endParameter(e, b), at(b));
            }
        }
        return false;
    }
    if (kf == Location::Edge && kt == Location::Point) {
        const auto ends = adjacency_.edgeVertices(a);
        if (ends[0] == b || ends[1] == b) {
            return add(-1, a, ends[0] == b ? 0.0f : 1.0f, at(b));
        }
        return false;
    }
    if (kf == Location::Point && kt == Location::Edge) {
        const auto ends = adjacency_.edgeVertices(b);
        if (ends[0] == a || ends[1] == a) {
            const float d0 = distance(to.position, at(ends[0]));
            const float d1 = distance(at(ends[0]), at(ends[1]));
            appendEntry(out, -1, b, std::min(d0 / d1, 1.0f), to.position);
            return true;
        }
        return false;
    }
    return false;
}

bool Tracer::candidates(const Location& here, const Spot& target, std::vector<int>& halfEdges,
                        CurveSegment& out) const
{
    halfEdges.clear();
    const Location& goal = target.location;
    if (here.kind == Location::Face) {
        if (goal.kind == Location::Face && goal.index == here.index) {
            return true;
        }
        for (int c = topology_.faceBegin(here.index); c < topology_.faceEnd(here.index); ++c) {
            halfEdges.push_back(c);
        }
        return false;
    }
    if (here.kind == Location::Edge) {
        const int e = here.index;
        const auto ends = adjacency_.edgeVertices(e);
        if (goal.kind == Location::Point && (ends[0] == goal.index || ends[1] == goal.index)) {
            // Along the edge to its end point: the face of the first corner of the edge.
            const auto corners = adjacency_.edgeCorners(e);
            const int face = corners.empty() ? -1 : adjacency_.cornerFace(corners.front());
            appendPathEntry(out, face, e, endParameter(e, goal.index), at(goal.index));
            return true;
        }
        for (const int c : adjacency_.edgeCorners(e)) {
            const int face = adjacency_.cornerFace(c);
            if (goal.kind == Location::Face && goal.index == face) {
                return true;
            }
            // Every other side of the face; the target among their start corners.
            const int m = topology_.faceSize(face);
            for (int r = 0; r + 1 < m; ++r) {
                const int k = cornerAfter(c, r + 1);
                halfEdges.push_back(k);
                if (r != 0 && goal.kind == Location::Point && goal.index == adjacency_.vertex(k)) {
                    const int edge = firstEdge(goal.index);
                    appendPathEntry(out, face, edge, endParameter(edge, goal.index), at(goal.index));
                    return true;
                }
            }
        }
        return false;
    }
    if (here.kind == Location::Point) {
        for (const int c : adjacency_.pointCorners(here.index)) {
            const int face = adjacency_.cornerFace(c);
            if (goal.kind == Location::Face && goal.index == face) {
                return true;
            }
            const int m = topology_.faceSize(face);
            if (goal.kind == Location::Point) {
                // The target among the other corners of the face: along a side
                // when it neighbours the point, else across the face. A
                // triangle mesh looks only at the next corner.
                const int looked = topology_.mode() == FaceMode::Triangles ? 1 : m - 1;
                for (int r = 0; r < looked; ++r) {
                    if (goal.index == adjacency_.vertex(cornerAfter(c, r + 1))) {
                        const int edge = firstEdge(goal.index);
                        const int across = r == 0 || r == m - 2 ? -1 : face;
                        appendPathEntry(out, across, edge, endParameter(edge, goal.index), at(goal.index));
                        return true;
                    }
                }
            }
            // The sides of the face that do not touch the point.
            for (int r = 1; r + 1 < m; ++r) {
                halfEdges.push_back(cornerAfter(c, r));
            }
        }
    }
    return false;
}

Approach Tracer::approach(const Vec3d& c, const Vec3d& direction, Vec3f a, Vec3f b) const
{
    Approach out;
    Vec3d u = wide(b) - wide(a);
    const double edgeLength = normalizeInPlace(u);
    const Vec3d w0 = c - wide(a);
    const Vec3d crossed = cross(u, direction);
    const double crossLength = length(crossed);
    const auto w = static_cast<float>(dot(u, direction));
    const auto den = static_cast<float>(-1.0f + w * w);
    if ((-1e-5 < crossLength && crossLength < 1e-5) || den == 0.0f) {
        out.parallel = true;
        out.s = dot(direction, wide(a) - c);
        out.t = 0.0;
        out.distance = length(cross(w0, direction));
        return out;
    }
    const auto along = static_cast<float>(dot(w0, direction));
    const auto acrossEdge = static_cast<float>(dot(u, w0));
    const float numS = along - w * acrossEdge;
    const float numT = w * along - acrossEdge;
    out.s = static_cast<double>(numS) / static_cast<double>(den);
    out.t = static_cast<double>(numT) / static_cast<double>(den) / edgeLength;
    return out;
}

// Corner fallback of the straight walk when no crossing lies ahead: the end
// points v of the candidate edges (in half edge order), with q their position
// along `direction` from c and g their distance to that line. `best` holds
// the result of the last pass and is replaced on every choice by v (s = q,
// t = 0/1 on its edge):
//  - after a choice with s > 0: v when q > 0 and g is smaller than the best g
//    (taken over only from choices with q > 0),
//  - at a point: v when q > s,
//  - on an edge or in a face: v when s < -0.01; before any choice, on the
//    edge of the crossing only v0 with t < 0.01 or v1 with t > 0.99; else
//    when q > s and, with L = |crossing - c| rounded to float,
//    (L < 1e-4 and g < 1e-4) or g < (L > 1e-5 ? 100 L : 1e-3).
void Tracer::cornerFallback(const Location& here, const Vec3d& c, const Vec3d& direction,
                            const std::vector<int>& halfEdges, Crossing& best) const
{
    const bool atPoint = here.kind != Location::Edge && here.kind != Location::Face;
    bool chosen = false;
    double bestGap = static_cast<double>(FLT_MAX);
    for (const int h : halfEdges) {
        const int e = adjacency_.cornerEdge(h);
        if (e < 0) {
            continue;
        }
        const auto ends = adjacency_.edgeVertices(e);
        for (const int corner : {h, adjacency_.next(h)}) {
            const int v = adjacency_.vertex(corner);
            const Vec3d r = wide(at(v)) - c;
            const double q = dot(r, direction);
            const Vec3d off{r.x - q * direction.x, r.y - q * direction.y, r.z - q * direction.z};
            const double gap = length(off);
            bool take = false;
            if (best.s > 0.0) {
                take = q > 0.0 && bestGap > gap;
            } else if (atPoint) {
                take = q > best.s;
            } else if (-0.01 > best.s) {
                take = true;
            } else if (!chosen && e == best.edge) {
                take = (0.01 > best.t && ends[0] == v) || (best.t > 0.99 && ends[1] == v);
            } else if (q > best.s) {
                const double reach = wide(static_cast<float>(length(best.point - c)));
                if (1e-4 > reach && 1e-4 > gap) {
                    take = true;
                } else {
                    take = (reach > 1e-5 ? reach * 100.0 : 1e-3) > gap;
                }
            }
            if (!take) {
                continue;
            }
            chosen = true;
            if (q > 0.0) {
                bestGap = gap;
            }
            best = {true, q, ends[0] == v ? 0.0 : 1.0, e, adjacency_.cornerFace(h), wide(at(v))};
        }
    }
}

bool Tracer::walk(const Spot& start, const Spot& target, CurveSegment& out, bool& backward) const
{
    backward = false;
    Vec3d c = wide(start.position);
    Location here = start.location;
    Location before{};
    const Vec3d goal = wide(target.position);
    std::vector<int> halfEdges;
    // Window widths are float products of the tolerance, then widened.
    const double a = absoluteTolerance;
    const double a3 = 3.0f * absoluteTolerance;
    const double a10 = 10.0f * absoluteTolerance;
    const double a100 = 100.0f * absoluteTolerance;
    for (int step = 0; step < maxWalkSteps; ++step) {
        if (candidates(here, target, halfEdges, out)) {
            return true;
        }
        if (c.x == goal.x && c.y == goal.y && c.z == goal.z) {
            backward = false;
            return false;
        }
        Vec3d d = goal - c;
        normalizeInPlace(d);
        // Projected direction of the last candidate examined (with an edge);
        // the corner fallback and the choice of an edge end measure along it.
        Vec3d lastDirection{0.0, 0.0, 1.0};
        auto pass = [&](double minLength, double lowFrom, double lowTo, double highFrom, double highTo) {
            Crossing best;
            for (const int h : halfEdges) {
                const int e = adjacency_.cornerEdge(h);
                if (e < 0) {
                    continue;
                }
                const int face = adjacency_.cornerFace(h);
                const Vec3d n = wide((*normals_)[static_cast<std::size_t>(face)]);
                const double dn = dot(d, n);
                Vec3d dp{d.x - dn * n.x, d.y - dn * n.y, d.z - dn * n.z};
                const double projectedLength = normalizeInPlace(dp);
                lastDirection = dp;
                if (!(minLength <= projectedLength)) {
                    continue;
                }
                const auto ends = adjacency_.edgeVertices(e);
                const Vec3d p0 = wide(at(ends[0]));
                const Vec3d p1 = wide(at(ends[1]));
                Approach hit = approach(c, dp, at(ends[0]), at(ends[1]));
                double t = hit.t;
                if (!hit.parallel) {
                    if (t > 0.98) {
                        const Approach reversed = approach(c, dp, at(ends[1]), at(ends[0]));
                        hit.s = reversed.s;
                        t = 1.0 - reversed.t;
                    }
                    if (lowFrom <= t && t <= lowTo) {
                        t = 0.0;
                    } else if (highFrom <= t && t <= highTo) {
                        t = 1.0;
                    } else if (!(0.0 <= t && t <= 1.0)) {
                        continue;
                    }
                } else {
                    // Distance rounded to float.
                    if (wide(static_cast<float>(hit.distance)) > 1e-6) {
                        continue;
                    }
                    // The end of the edge that lies further along the direction,
                    // the second one on ties; s is its position along it.
                    const double q0 = dot(p0 - c, dp);
                    const double q1 = dot(p1 - c, dp);
                    const bool first = q0 > q1;
                    t = first ? 0.0 : 1.0;
                    hit.s = first ? q0 : q1;
                }
                if (hit.s > best.s) {
                    const double r = 1.0 - t;
                    best = {true, hit.s, t, e, face, {r * p0.x + t * p1.x, r * p0.y + t * p1.y, r * p0.z + t * p1.z}};
                }
            }
            return best;
        };
        Crossing best = pass(0.3, -a3, a, 1.0 - a, 1.0 + a3);
        if (!best.valid || best.s < 0.0) {
            best = pass(0.05, -a100, a10, 1.0 - a10, 1.0 + a100);
        }
        if (!best.valid || best.s < 0.0) {
            cornerFallback(here, c, lastDirection, halfEdges, best);
            if (!best.valid) {
                backward = false;
                return false;
            }
        }
        if (best.s < 0.0) {
            const float limit = 5.0f * absoluteTolerance;
            if (here.kind == Location::Edge && !out.empty() && out.edges.back() != -1 &&
                (out.t.back() < limit || out.t.back() > 1.0f - limit)) {
                // The last entry moves onto its corner point; no new entry.
                const int lastEdge = out.edges.back();
                const auto lastEnds = adjacency_.edgeVertices(lastEdge);
                const bool first = out.t.back() < limit;
                out.t.back() = first ? 0.0f : 1.0f;
                const int v = first ? lastEnds[0] : lastEnds[1];
                out.points.back() = at(v);
                // The location before the current one stays as it was.
                here = {Location::Point, v};
                c = wide(at(v));
                if (here == target.location) {
                    return true;
                }
                continue;
            }
            backward = true;
            if (here.kind == Location::Edge) {
                // The end of the current edge whose position along the last
                // projected direction lies nearer to that of the target (the
                // second one on ties).
                const auto ends = adjacency_.edgeVertices(here.index);
                const double q0 = dot(wide(at(ends[0])) - c, lastDirection);
                const double q1 = dot(wide(at(ends[1])) - c, lastDirection);
                const double qz = dot(goal - c, lastDirection);
                const bool first = std::fabs(q1 - qz) > std::fabs(q0 - qz);
                best.edge = here.index;
                best.t = first ? 0.0 : 1.0;
                best.s = first ? q0 : q1;
                best.point = wide(at(first ? ends[0] : ends[1]));
            }
        }
        // Enter the crossing: snapping to an end point changes t only, the
        // crossing point itself is kept (in double, also as the next c).
        const auto ends = adjacency_.edgeVertices(best.edge);
        const Vec3d d0 = best.point - wide(at(ends[0]));
        const Vec3d d1 = best.point - wide(at(ends[1]));
        const double snap2 = wide(tolerance2_);
        const bool near0 = dot(d0, d0) < snap2;
        const bool near1 = dot(d1, d1) < snap2;
        double t = best.t;
        if (near0 && near1) {
            t = t > 0.9 ? 1.0 : 0.0;
        } else if (near0) {
            t = 0.0;
        } else if (near1) {
            t = 1.0;
        }
        Location next{Location::Edge, best.edge};
        if (t == 0.0) {
            next = {Location::Point, ends[0]};
        } else if (t == 1.0) {
            next = {Location::Point, ends[1]};
        }
        if (step > 0 && next == before) {
            out.clear(); // back and forth
            backward = false;
            return false;
        }
        const Vec3f hitPoint{static_cast<float>(best.point.x), static_cast<float>(best.point.y),
                             static_cast<float>(best.point.z)};
        appendPathEntry(out, best.face, best.edge, static_cast<float>(t), hitPoint);
        c = best.point;
        before = here;
        here = next;
        if (here == target.location) {
            return true;
        }
    }
    out.clear();
    backward = false;
    return false;
}

bool Tracer::edgePath(const Spot& start, const Spot& target, CurveSegment& out) const
{
    out.clear();
    const int pointCount = mesh_.pointCount();
    const Vec3f s = start.position;
    float bound = 20.0f * distance(s, target.position);
    std::vector<float> cost(static_cast<std::size_t>(pointCount), bound);
    std::vector<int> previous(static_cast<std::size_t>(pointCount), -1);
    std::deque<int> queue;
    auto pointsOf = [&](const Location& l) {
        std::vector<int> list;
        if (l.kind == Location::Point) {
            list.push_back(l.index);
        } else if (l.kind == Location::Edge) {
            const auto ends = adjacency_.edgeVertices(l.index);
            list.assign(ends.begin(), ends.end());
        } else if (l.kind == Location::Face) {
            for (int c = 3 * l.index; c < 3 * l.index + 3; ++c) {
                list.push_back(topology_.vertex(c));
            }
        }
        return list;
    };
    for (const int v : pointsOf(target.location)) {
        cost[static_cast<std::size_t>(v)] = target.location.kind == Location::Point ? 0.0f : distance(at(v), s);
        previous[static_cast<std::size_t>(v)] = -2;
        queue.push_back(v);
    }
    const std::vector<int> starts = pointsOf(start.location);
    for (const int v : starts) {
        previous[static_cast<std::size_t>(v)] = -3;
    }
    std::vector<int> order;
    std::vector<float> keys;
    while (!queue.empty()) {
        const int u = queue.front();
        queue.pop_front();
        const auto neighbours = adjacency_.pointNeighbours(u);
        keys.clear();
        for (const int w : neighbours) {
            const Vec3f d = at(w) - s;
            keys.push_back(dot(d, d));
        }
        order.resize(neighbours.size());
        std::iota(order.begin(), order.end(), 0);
        util::Introsort::sort(order.begin(), order.end(), [&](int x, int y) {
            return keys[static_cast<std::size_t>(x)] < keys[static_cast<std::size_t>(y)];
        });
        for (const int i : order) {
            const int w = neighbours[static_cast<std::size_t>(i)];
            const float c = cost[static_cast<std::size_t>(u)] + edgeLength(at(u), at(w));
            if (previous[static_cast<std::size_t>(w)] == -3) {
                bound = c;
            } else if (c < cost[static_cast<std::size_t>(w)] && c < bound) {
                cost[static_cast<std::size_t>(w)] = c;
                previous[static_cast<std::size_t>(w)] = u;
                if (c + distance(at(w), s) < bound) {
                    queue.push_back(w);
                }
            }
        }
    }
    // Every start point takes its cheapest reached neighbour.
    for (const int g : starts) {
        for (const int w : adjacency_.pointNeighbours(g)) {
            const int pw = previous[static_cast<std::size_t>(w)];
            if (pw == -1 || pw == -3) {
                continue;
            }
            const float c = cost[static_cast<std::size_t>(w)] + edgeLength(at(g), at(w));
            if (c < cost[static_cast<std::size_t>(g)]) {
                cost[static_cast<std::size_t>(g)] = c;
                previous[static_cast<std::size_t>(g)] = w;
            }
        }
    }
    int chosen = -1;
    if (start.location.kind == Location::Point) {
        chosen = start.location.index;
    } else {
        float best = FLT_MAX;
        for (const int g : starts) {
            if (previous[static_cast<std::size_t>(g)] < 0) {
                continue;
            }
            const float value = distance(at(g), s) + cost[static_cast<std::size_t>(g)];
            if (value < best) {
                best = value;
                chosen = g;
            }
        }
    }
    if (chosen < 0) {
        return false;
    }
    auto push = [&](int v) {
        const int e = firstEdge(v);
        out.edges.push_back(e);
        out.t.push_back(e >= 0 ? endParameter(e, v) : 0.0f);
        out.points.push_back(at(v));
        out.faces.push_back(-1);
    };
    if (start.location.kind != Location::Point) {
        push(chosen);
    }
    int current = chosen;
    while (true) {
        const int next = previous[static_cast<std::size_t>(current)];
        if (next == -2) {
            break;
        }
        if (next < 0 || out.size() > pointCount) {
            out.clear();
            return false;
        }
        push(next);
        current = next;
    }
    return true;
}

float Tracer::pathLength(const Spot& start, const CurveSegment& path, const Spot& target) const
{
    if (path.empty()) {
        return distance(start.position, target.position);
    }
    float sum = distance(start.position, path.points.front());
    for (std::size_t k = 0; k + 1 < path.points.size(); ++k) {
        sum += distance(path.points[k], path.points[k + 1]);
    }
    if (target.location.kind != Location::Point) {
        sum += distance(target.position, path.points.back());
    }
    return sum;
}

// Mean length of the edges at the start and then at the target, each only
// when it is a point (edge lengths as the mesh measures keep them, float sum
// in point -> edge order); 0 without edges.
float Tracer::chordLimit(const Spot& start, const Spot& target) const
{
    float sum = 0.0f;
    int count = 0;
    for (const Spot* spot : {&start, &target}) {
        if (spot->location.kind != Location::Point) {
            continue;
        }
        for (const int e : adjacency_.pointEdges(spot->location.index)) {
            const auto ends = adjacency_.edgeVertices(e);
            sum += edgeLength(at(ends[0]), at(ends[1]));
            ++count;
        }
    }
    return count != 0 ? sum / static_cast<float>(count) : sum;
}

// Moves the points of a path onto the chord start -> target, in order. With
// `giveUp` every point moves and the first one further than the limit stops
// it (false; the points before stay moved); otherwise only points nearer than
// the limit move.
bool Tracer::layOntoChord(const Spot& start, const Spot& target, CurveSegment& path, bool giveUp) const
{
    const float limit = chordLimit(start, target);
    for (Vec3f& p : path.points) {
        const SegmentClosest closest = closestOnSegment(p, start.position, target.position);
        const float gap = std::sqrt(closest.distanceSq);
        if (giveUp) {
            if (limit < gap) {
                return false;
            }
            p = closest.point;
        } else if (gap < limit) {
            p = closest.point;
        }
    }
    return true;
}

std::vector<CurveSegment> Tracer::project(const LocatedPolyline& curve) const
{
    std::vector<CurveSegment> out;
    const std::size_t n = std::min(curve.points.size(), curve.locations.size());
    if (n == 0) {
        return out;
    }
    const bool closed = n > 1 && curve.points.front() == curve.points[n - 1] && curve.locations.front() == curve.locations[n - 1];
    CurveSegment piece;
    piece.points.reserve(50);
    Spot previous{curve.points[0], curve.locations[0]};
    std::size_t previousIndex = 0;
    // Codes of the entries appended last: -1, the last one reaches point `index`.
    auto addCodes = [&](std::size_t index) {
        if (!options_.sourceCodes) {
            return;
        }
        while (piece.sources.size() < piece.edges.size()) {
            piece.sources.push_back(piece.sources.size() + 1 == piece.edges.size() ? static_cast<int>(index) : -1);
        }
    };
    locate(previous);
    auto writeFirst = [&](const Spot& spot) {
        if (spot.location.kind == Location::Face) {
            appendEntry(piece, -1, -1, 0.0f, spot.position);
        } else if (spot.location.kind == Location::Point) {
            const int e = firstEdge(spot.location.index);
            appendEntry(piece, -1, e, e >= 0 ? endParameter(e, spot.location.index) : 0.0f, at(spot.location.index));
        } else if (spot.location.kind == Location::Edge) {
            const auto ends = adjacency_.edgeVertices(spot.location.index);
            const SegmentPoint q = nearestOnSegment(spot.position, at(ends[0]), at(ends[1]));
            appendEntry(piece, -1, spot.location.index, q.t, q.point);
        }
    };
    writeFirst(previous);
    addCodes(0);
    for (std::size_t i = 1; i < n; ++i) {
        Spot next{curve.points[i], curve.locations[i]};
        if (previous.location.kind == Location::Point && next.location == previous.location) {
            continue;
        }
        locate(next);
        if (connect(previous, next, piece)) {
            addCodes(i);
            previous = next;
            previousIndex = i;
            continue;
        }
        CurveSegment straight;
        bool backward = false;
        const bool walked = walk(previous, next, straight, backward);
        const bool kept = walked && (!options_.chordCheck || layOntoChord(previous, next, straight, true));
        const CurveSegment* path = nullptr;
        CurveSegment edges;
        if (kept && !backward) {
            path = &straight;
        } else {
            if (straight.empty()) {
                backward = false;
            }
            if (!edgePath(previous, next, edges)) {
                // Neither reaches the point: the piece ends here.
                if (piece.size() >= 2) {
                    out.push_back(std::move(piece));
                }
                piece = CurveSegment{};
                previous = next;
                previousIndex = i;
                continue;
            }
            if (options_.chordCheck) {
                (void)layOntoChord(previous, next, edges, false);
            }
            if (!backward) {
                path = &edges;
            } else {
                const float lb = pathLength(previous, straight, next);
                const float lw = pathLength(previous, edges, next);
                path = lb >= lw ? &edges : &straight;
            }
        }
        if (!path->empty()) {
            const bool wasEmpty = piece.empty();
            piece.edges.insert(piece.edges.end(), path->edges.begin(), path->edges.end());
            piece.t.insert(piece.t.end(), path->t.begin(), path->t.end());
            piece.points.insert(piece.points.end(), path->points.begin(), path->points.end());
            piece.faces.insert(piece.faces.end(), path->faces.begin() + (wasEmpty ? 1 : 0), path->faces.end());
            addCodes(i);
        }
        previous = next;
        previousIndex = i;
    }
    if (piece.empty()) {
        return out;
    }
    if (previous.location.kind == Location::Face && !piece.faces.empty()) {
        const bool same = piece.faces.back() == previous.location.index && piece.points.back() == previous.position;
        if (!same) {
            appendEntry(piece, previous.location.index, -1, 0.0f, previous.position);
            addCodes(previousIndex);
        }
    }
    if (closed) {
        const int e0 = piece.edges.front();
        const int e1 = piece.edges.back();
        const float t0 = piece.t.front();
        const float t1 = piece.t.back();
        if (!(e0 == e1 && t0 == t1)) {
            if (e0 == e1 && static_cast<double>(std::fabs(t0 - t1)) < 1e-6) {
                piece.t.back() = t0;
                piece.points.back() = piece.points.front();
            } else if (e0 >= 0 && e1 >= 0) {
                const int f = commonFace(e0, e1);
                if (f >= 0) {
                    appendEntry(piece, f, e0, t0, piece.points.front());
                    if (options_.sourceCodes) {
                        piece.sources.push_back(piece.sources.front());
                    }
                }
            }
        }
    }
    out.push_back(std::move(piece));
    return out;
}

} // namespace

std::vector<CurveSegment> CurveProjector::project(const Mesh& mesh, const LocatedPolyline& curve)
{
    return project(mesh, curve, Options{});
}

std::vector<CurveSegment> CurveProjector::project(const Mesh& mesh, const LocatedPolyline& curve, const Options& options)
{
    return Tracer(mesh, options).project(curve);
}

std::vector<std::vector<CurveSegment>> CurveProjector::projectAll(const Mesh& mesh, std::span<const LocatedPolyline> curves,
                                                                  const Options& options)
{
    const Tracer tracer(mesh, options);
    std::vector<std::vector<CurveSegment>> out;
    out.reserve(curves.size());
    for (const LocatedPolyline& curve : curves) {
        out.push_back(tracer.project(curve));
    }
    return out;
}

std::vector<int> CurveProjector::edgePath(const Mesh& mesh, int from, int to)
{
    const Tracer tracer(mesh, Options{}, false);
    CurveSegment path;
    std::vector<int> out;
    const Spot start{mesh.points()[static_cast<std::size_t>(from)], {Location::Point, from}};
    const Spot target{mesh.points()[static_cast<std::size_t>(to)], {Location::Point, to}};
    if (!tracer.edgePath(start, target, path)) {
        return out;
    }
    const Adjacency& adjacency = mesh.adjacency();
    for (int i = 0; i < path.size(); ++i) {
        const int e = path.edges[static_cast<std::size_t>(i)];
        if (e < 0) {
            continue;
        }
        const auto ends = adjacency.edgeVertices(e);
        out.push_back(path.t[static_cast<std::size_t>(i)] == 0.0f ? ends[0] : ends[1]);
    }
    return out;
}

} // namespace areko::remesher::curves
