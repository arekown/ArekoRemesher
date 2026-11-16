#include "symmetry_detector.h"

#include "mesh/measures.h"
#include "mesh/normals.h"
#include "mesh/vertex_fan.h"
#include "util/index_queue.h"
#include "util/introsort.h"

#include <array>
#include <bit>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr float scoreLimit = std::bit_cast<float>(0x3f7d70a4u);      // 0.99
constexpr float strictDirection = std::bit_cast<float>(0x3a83126fu); // 0.001
constexpr float directionWeight = std::bit_cast<float>(0x3f333333u); // 0.7, also the share of paired faces
constexpr float positionWeight = std::bit_cast<float>(0x3e99999au);  // 0.3
constexpr float snapLimit = std::bit_cast<float>(0x3f7eb852u);       // 0.995
constexpr float strictPairDivisor = 200.0f;
constexpr int changeLimit = 10;
constexpr std::size_t longLine = 500; // from this many entries the limit is the length / 50
constexpr int longLineDivisor = 50;
constexpr float narrowDivisor = 500.0f;      // narrow: below the bounding box diagonal / 500
constexpr double leavingShare = 0.05;        // share along the axis from which a line leaves the plane
constexpr float bendLimit = degF(10.0);      // largest turn of a line in the tangent plane
static_assert(std::bit_cast<std::uint32_t>(bendLimit) == 0x3e32b8c2u);

// Start edge limits of one round: score above `score` and at most
// `scoreMax`, direction share along the axis at most `direction`, distance
// from the middle at most `position` (relative to half the width).
struct Round {
    float score;
    float scoreMax;
    float direction;
    float position;
};
// Further rounds with the same limits would find no new start edges.
constexpr std::array<Round, 2> rounds{{{scoreLimit, 1.0f, strictDirection, 1.0f}, {scoreLimit, 1.0f, 1.0f, 1.0f}}};

enum class Kind { Point, Edge, Face };

// Line entries come in pairs: (edge, point) along an edge to its end point,
// (face, edge) through a face over an edge, (face, point) through a face to
// its opposite corner.
struct Entry {
    Kind kind;
    int index;
};

// Straightness of one followed line: the start edge is narrow along an axis
// when its extent there is below the limit; the check waits while one axis
// stays narrow. `frame` is the normal of the tangent plane, zero until known.
struct Straightness {
    std::array<bool, 3> narrow{};
    int narrowCount = 0;
    Vec3f frame{};
};

// Pairs found from one start edge.
struct Pairing {
    std::vector<int> faces;
    std::vector<int> edges;
    std::vector<int> points;
    BitArray side; // faces of the searched side
    int paired = 0;
    int lineLength = 0;
};

bool contradicts(const std::vector<int>& mirror, int p, int q)
{
    const int mp = mirror[static_cast<std::size_t>(p)];
    const int mq = mirror[static_cast<std::size_t>(q)];
    return (mp != -1 && mp != q) || (mq != -1 && mq != p);
}

void pair(std::vector<int>& mirror, int u, int w)
{
    mirror[static_cast<std::size_t>(u)] = w;
    mirror[static_cast<std::size_t>(w)] = u;
}

class Matcher {
public:
    Matcher(const Mesh& mesh, const SymmetryCache* previous)
        : mesh_(mesh), adjacency_(mesh.adjacency()), points_(mesh.points()), faceCount_(mesh.faceCount()),
          narrowLimit_(diagonal(mesh.bounds()) / narrowDivisor),
          onLine_(static_cast<std::size_t>(adjacency_.edgeCount())),
          crossed_(static_cast<std::size_t>(adjacency_.edgeCount()))
    {
        if (previous != nullptr) {
            lock_ = previous->lockedPlane();
        }
    }

    std::optional<SymmetryCache> run()
    {
        if (mesh_.topology().mode() != FaceMode::Triangles || faceCount_ == 0) {
            return std::nullopt;
        }
        Pairing best;
        int bestPaired = -1;
        int bestAxis = -1;
        bool complete = false;
        for (const Round& round : rounds) {
            const std::array<std::vector<int>, 3> lists = startEdges(round);
            // Axes by list size ascending (stable), taken from the back.
            std::array<int, 3> order{0, 1, 2};
            util::Introsort::sort(order.begin(), order.end(), [&](int a, int b) {
                return lists[static_cast<std::size_t>(a)].size() < lists[static_cast<std::size_t>(b)].size();
            });
            for (int i = 2; i >= 0 && !complete; --i) {
                const int axis = order[static_cast<std::size_t>(i)];
                for (const int e : lists[static_cast<std::size_t>(axis)]) {
                    if (onLine_.test(static_cast<std::size_t>(e))) {
                        continue;
                    }
                    Pairing candidate = pairFaces(e, axis);
                    if (candidate.paired <= bestPaired || candidate.lineLength < 1) {
                        continue;
                    }
                    best = std::move(candidate);
                    bestAxis = axis;
                    if (best.paired == faceCount_ || best.paired > faceCount_ - best.lineLength / 2) {
                        complete = true;
                        break;
                    }
                    bestPaired = best.paired;
                }
            }
            if (complete) {
                break;
            }
        }
        if (best.paired <= truncToInt(static_cast<float>(faceCount_) * directionWeight)) {
            return std::nullopt;
        }
        return take(std::move(best), bestAxis);
    }

private:
    [[nodiscard]] Vec3f at(int v) const { return points_[static_cast<std::size_t>(v)]; }
    [[nodiscard]] int face(int c) const { return adjacency_.cornerFace(c); }

    [[nodiscard]] const SplitNormals& splitNormals()
    {
        if (!splitNormals_) {
            splitNormals_.emplace(mesh_);
        }
        return *splitNormals_;
    }

    // Where a line entry lies: the middle of an edge.
    [[nodiscard]] Vec3f edgePoint(int e) const
    {
        const auto ends = adjacency_.edgeVertices(e);
        return (at(ends[0]) + at(ends[1])) * 0.5f;
    }

    // Centre of triangle f: ((0 + p0) + p1) + p2, divided by 3.
    [[nodiscard]] Vec3f faceCentre(int f) const
    {
        const Topology& topology = mesh_.topology();
        Vec3f sum{};
        for (int k = 0; k < 3; ++k) {
            sum = sum + at(topology.vertex(3 * f + k));
        }
        return {sum.x / 3.0f, sum.y / 3.0f, sum.z / 3.0f};
    }

    // Straightness check after a step (at least four entries): with an axis
    // still narrow the check waits; each edge step that is wide along a
    // narrow axis ends the wait for that axis. Then the points A (end point
    // of the entry edge or face centre), B and C (points or edge middles) of
    // the last two steps give u = B - A and w = C - B. When one of them leaves
    // the axis plane (share at least 0.05), both are projected into the
    // tangent plane and the line bends when they turn by more than the limit.
    // The tangent plane comes from the entry edge, else from B's or C's edge
    // while still unknown; without it nothing is checked.
    bool bends(const std::vector<Entry>& line, int axis, Straightness& s)
    {
        const std::size_t n = line.size();
        const Entry& a = line[n - 4];
        const Entry& b = line[n - 3];
        const Entry& step = line[n - 2];
        const Entry& c = line[n - 1];
        if (step.kind == Kind::Edge && s.narrowCount > 0) {
            const auto ends = adjacency_.edgeVertices(step.index);
            const Vec3f d = at(ends[0]) - at(ends[1]);
            for (int k = 0; k < 3; ++k) {
                if (s.narrow[static_cast<std::size_t>(k)] && std::fabs(axisValue(d, k)) >= narrowLimit_) {
                    s.narrow[static_cast<std::size_t>(k)] = false;
                    --s.narrowCount;
                }
            }
            if (s.narrowCount > 0) {
                return false;
            }
        }
        Vec3f pa;
        Vec3f pb;
        if (a.kind == Kind::Edge) {
            const auto ends = adjacency_.edgeVertices(a.index);
            pa = at(b.index == ends[0] ? ends[1] : ends[0]);
            pb = at(b.index);
            s.frame = splitNormals().edgeNormal(adjacency_, a.index);
        } else {
            pa = faceCentre(a.index);
            if (b.kind == Kind::Point) {
                pb = at(b.index);
            } else {
                pb = edgePoint(b.index);
                if (s.frame == Vec3f{}) {
                    s.frame = splitNormals().edgeNormal(adjacency_, b.index);
                }
            }
        }
        Vec3f pc;
        if (c.kind == Kind::Edge) {
            pc = edgePoint(c.index);
            if (s.frame == Vec3f{}) {
                s.frame = splitNormals().edgeNormal(adjacency_, c.index);
            }
        } else {
            pc = at(c.index);
        }
        if (s.frame == Vec3f{} || axis < 0) {
            return false;
        }
        Vec3f u = pb - pa;
        Vec3f w = pc - pb;
        const float shareU = std::fabs(axisValue(u, axis) / length(u));
        const float shareW = std::fabs(axisValue(w, axis) / length(w));
        if (leavingShare > wide(shareU) && leavingShare > wide(shareW)) {
            return false;
        }
        u = projectOntoPlane(u, s.frame);
        w = projectOntoPlane(w, s.frame);
        return lineAngle(u, w) > bendLimit;
    }

    // Edges near the middle plane of each axis, ascending; an edge joins the
    // first axis it fits.
    [[nodiscard]] std::array<std::vector<int>, 3> startEdges(const Round& round) const
    {
        const Bounds& box = mesh_.bounds();
        const int only = lock_ && lock_->axis >= 0 && lock_->axis <= 2 ? lock_->axis : -1;
        std::array<std::vector<int>, 3> lists;
        for (int e = 0; e < adjacency_.edgeCount(); ++e) {
            if (onLine_.test(static_cast<std::size_t>(e))) {
                continue;
            }
            const auto ends = adjacency_.edgeVertices(e);
            const Vec3f p0 = at(ends[0]);
            const Vec3f p1 = at(ends[1]);
            Vec3f d = p1 - p0;
            normalizeInPlace(d);
            for (int a = 0; a < 3; ++a) {
                if (only >= 0 && a != only) {
                    continue;
                }
                const float share = std::fabs(axisValue(d, a));
                if (share > round.direction) {
                    continue;
                }
                const float lo = axisValue(box.min, a);
                const float hi = axisValue(box.max, a);
                const float width = hi - lo;
                if (width == 0.0f) {
                    continue;
                }
                const float middle = (lo + hi) * 0.5f;
                const float r0 = std::fabs(axisValue(p0, a) - middle);
                const float r1 = std::fabs(axisValue(p1, a) - middle);
                const float r = (r0 > r1 ? r0 : r1) / (width * 0.5f);
                if (r > round.position) {
                    continue;
                }
                const float score = (1.0f - share) * directionWeight + (1.0f - r) * positionWeight;
                if (!(score > round.score) || round.scoreMax < score) {
                    continue;
                }
                lists[static_cast<std::size_t>(a)].push_back(e);
                break;
            }
        }
        return lists;
    }

    // One step along edge `next` from `v`; false when the edge is missing or
    // already on a line.
    bool stepAlong(int next, int& v, std::vector<Entry>& line)
    {
        if (next < 0 || onLine_.test(static_cast<std::size_t>(next))) {
            return false;
        }
        onLine_.set(static_cast<std::size_t>(next));
        v = adjacency_.otherEnd(next, v);
        line.push_back({Kind::Edge, next});
        line.push_back({Kind::Point, v});
        return true;
    }

    // Into the face of corner k0 and over its edge opposite k0; false when
    // that edge is missing or already crossed.
    bool crossFrom(int k0, std::vector<Entry>& line, int& inFace, int& inEdge)
    {
        const int f = face(k0);
        const int x = adjacency_.cornerEdge(adjacency_.next(k0));
        if (x < 0 || crossed_.test(static_cast<std::size_t>(x))) {
            return false;
        }
        crossed_.set(static_cast<std::size_t>(x));
        line.push_back({Kind::Face, f});
        line.push_back({Kind::Edge, x});
        inFace = f;
        inEdge = x;
        return true;
    }

    // Follows the line of edge e from its first (or second) end point. Faces
    // before the middle of each passed fan go to `sideA`, faces after it to
    // `sideB`. A line that bends away from the plane of `axis` is emptied and
    // ends. Returns the number of changes between edges and faces.
    int follow(int e, int axis, bool fromFirst, std::vector<Entry>& line, BitArray& sideA, BitArray& sideB)
    {
        const auto ends = adjacency_.edgeVertices(e);
        int v = fromFirst ? ends[0] : ends[1];
        int changes = 0;
        Straightness straightness;
        const Vec3f extent = at(ends[0]) - at(ends[1]);
        for (int k = 0; k < 3; ++k) {
            const bool narrow = std::fabs(axisValue(extent, k)) < narrowLimit_;
            straightness.narrow[static_cast<std::size_t>(k)] = narrow;
            straightness.narrowCount += narrow ? 1 : 0;
        }
        line = {{Kind::Edge, e}, {Kind::Point, v}};
        onLine_.set(static_cast<std::size_t>(e));
        int inEdge = e;  // edge of the last step, -1 at an opposite corner
        int inFace = -1; // face of the last step, -1 at a point
        int corner = -1; // opposite corner reached in `inFace`
        while (true) {
            if (inFace < 0) {
                const VertexFan fan = vertexFan(adjacency_, v, inEdge);
                if (fan.corners.empty() || fan.edges.empty()) {
                    return changes;
                }
                const std::size_t m = fan.corners.size();
                if (fan.edges.front() != fan.edges.back()) {
                    if (line.size() == 2) {
                        for (const int k : fan.corners) {
                            sideA.setBounded(static_cast<std::size_t>(face(k)));
                        }
                    }
                    return changes;
                }
                if (m % 2 == 0) {
                    const int next = fan.edges[fan.edges.size() / 2];
                    if (!stepAlong(next, v, line)) {
                        return changes;
                    }
                    inEdge = next;
                } else {
                    ++changes;
                    if (!crossFrom(fan.corners[m / 2], line, inFace, inEdge)) {
                        return changes;
                    }
                }
                const std::size_t h = m / 2;
                for (std::size_t i = 0; i < h; ++i) {
                    sideA.setBounded(static_cast<std::size_t>(face(fan.corners[i])));
                }
                for (std::size_t i = m - h; i < m; ++i) {
                    sideB.setBounded(static_cast<std::size_t>(face(fan.corners[i])));
                }
            } else if (inEdge < 0) {
                const VertexFan fan = vertexFanFrom(adjacency_, corner);
                if (fan.edges.empty() || fan.edges.front() != fan.edges.back()) {
                    return changes;
                }
                const std::size_t m = fan.corners.size();
                if (m % 2 == 0) {
                    if (!crossFrom(fan.corners[m / 2], line, inFace, inEdge)) {
                        return changes;
                    }
                } else {
                    const int next = fan.edges[fan.edges.size() / 2];
                    if (!stepAlong(next, v, line)) {
                        return changes;
                    }
                    inFace = -1;
                    inEdge = next;
                    ++changes;
                }
            } else {
                // Over the crossed edge into the next face, to its corner opposite the edge.
                const auto corners = adjacency_.edgeCorners(inEdge);
                if (corners.size() != 2) {
                    return changes;
                }
                int other = -1;
                if (face(corners[0]) == inFace) {
                    other = corners[1];
                } else if (face(corners[1]) == inFace) {
                    other = corners[0];
                } else {
                    return changes;
                }
                inFace = face(other);
                corner = adjacency_.prev(other);
                v = adjacency_.vertex(corner);
                line.push_back({Kind::Face, inFace});
                line.push_back({Kind::Point, v});
                inEdge = -1;
            }
            if (line.size() >= 4 && bends(line, axis, straightness)) {
                line.clear();
                return changes;
            }
            const Entry& last = line[line.size() - 2];
            if (last.kind == Kind::Edge && last.index == e) {
                line.resize(line.size() - 2);
                return changes;
            }
        }
    }

    // The line through edge e in both directions: the first direction
    // reversed (without its start edge), then the second. `side` receives the
    // faces of one side. Returns the number of changes.
    int matchLine(int e, int axis, std::vector<Entry>& line, BitArray& side)
    {
        const auto edges = static_cast<std::size_t>(adjacency_.edgeCount());
        side = BitArray(edges);
        BitArray other(edges);
        std::vector<Entry> first;
        std::vector<Entry> second;
        const int c1 = follow(e, axis, true, first, side, other);
        const int c2 = follow(e, axis, false, second, other, side);
        line.clear();
        if (first.empty() && second.empty()) {
            return 0;
        }
        if (first.size() > 2) {
            line.assign(first.rbegin() + 1, first.rend() - 1);
        }
        line.insert(line.end(), second.begin(), second.end());
        return c1 + c2 - 1;
    }

    // Line of start edge e, then the breadth first pairing of faces from the
    // faces of the searched side at the line.
    Pairing pairFaces(int e, int axis)
    {
        std::vector<Entry> line;
        BitArray side;
        const int changes = matchLine(e, axis, line, side);
        const int limit = line.size() >= longLine ? static_cast<int>(line.size()) / longLineDivisor : changeLimit;
        Pairing r;
        if (changes > limit) {
            return r;
        }
        const auto faces = static_cast<std::size_t>(faceCount_);
        r.lineLength = static_cast<int>(line.size());
        r.faces.assign(faces, -1);
        r.edges.assign(static_cast<std::size_t>(adjacency_.edgeCount()), -1);
        r.points.assign(points_.size(), -1);
        r.side = BitArray(faces);
        util::IndexQueue queue(faces);
        std::vector<int> entry(faces, -1);

        // Faces on the line mirror onto themselves, edges along it too, and
        // a crossed edge swaps its end points.
        for (std::size_t i = 0; i + 1 < line.size(); i += 2) {
            const Entry& a = line[i];
            const Entry& b = line[i + 1];
            if (a.kind == Kind::Face) {
                if (r.faces[static_cast<std::size_t>(a.index)] == -1) {
                    ++r.paired;
                    r.faces[static_cast<std::size_t>(a.index)] = a.index;
                }
            } else if (a.kind == Kind::Edge) {
                r.edges[static_cast<std::size_t>(a.index)] = a.index;
                const auto ends = adjacency_.edgeVertices(a.index);
                r.points[static_cast<std::size_t>(ends[0])] = ends[0];
                r.points[static_cast<std::size_t>(ends[1])] = ends[1];
                for (const int k : adjacency_.edgeCorners(a.index)) {
                    const int g = face(k);
                    if (side.testBounded(static_cast<std::size_t>(g)) && !queue.contains(g)) {
                        entry[static_cast<std::size_t>(g)] = k;
                        queue.pushBack(g);
                    }
                }
            }
            if (b.kind == Kind::Edge) {
                r.edges[static_cast<std::size_t>(b.index)] = b.index;
                const auto ends = adjacency_.edgeVertices(b.index);
                pair(r.points, ends[0], ends[1]);
            }
        }

        std::vector<std::array<int, 2>> edgePairs;
        std::vector<std::array<int, 2>> pointPairs;
        std::vector<std::array<int, 2>> nextFaces;
        while (!queue.empty()) {
            const int f = queue.pop();
            if (r.faces[static_cast<std::size_t>(f)] >= 0) {
                continue;
            }
            const int k = entry[static_cast<std::size_t>(f)];
            const int x = k >= 0 ? adjacency_.cornerEdge(k) : -1;
            if (x < 0) {
                continue;
            }
            // Mirror face g across the mirror edge of the entry edge.
            int mirrorEdge = x;
            int g = -1;
            if (r.edges[static_cast<std::size_t>(x)] == x) {
                g = adjacency_.otherFace(x, f);
            } else {
                const int h = adjacency_.otherFace(x, f);
                mirrorEdge = r.edges[static_cast<std::size_t>(x)];
                const int mirrorOfH = h >= 0 ? r.faces[static_cast<std::size_t>(h)] : -1;
                g = mirrorEdge >= 0 && mirrorOfH >= 0 ? adjacency_.otherFace(mirrorEdge, mirrorOfH) : -1;
            }
            if (g < 0 || r.faces[static_cast<std::size_t>(g)] >= 0) {
                continue;
            }
            int kg = -1;
            for (const int c : adjacency_.edgeCorners(mirrorEdge)) {
                if (face(c) == g) {
                    kg = c;
                    break;
                }
            }
            if (kg < 0) {
                continue;
            }
            // Corners of f forwards against corners of g backwards.
            edgePairs.clear();
            pointPairs.clear();
            nextFaces.clear();
            bool fits = true;
            int a = adjacency_.next(k);
            int b = adjacency_.prev(kg);
            while (a != k) {
                const int ea = adjacency_.cornerEdge(a);
                if (ea >= 0) {
                    const int eb = adjacency_.cornerEdge(b);
                    const int known = r.edges[static_cast<std::size_t>(ea)];
                    if (known == -1) {
                        edgePairs.push_back({ea, eb});
                        const int pa = adjacency_.vertex(a);
                        const int pb = adjacency_.vertex(b);
                        const int pna = adjacency_.vertex(adjacency_.next(a));
                        const int pnb = adjacency_.vertex(adjacency_.next(b));
                        pointPairs.push_back({pa, pnb});
                        pointPairs.push_back({pna, pb});
                        if (contradicts(r.points, pa, pnb) || contradicts(r.points, pna, pb)) {
                            fits = false;
                            break;
                        }
                        const int across = adjacency_.otherFace(ea, f);
                        if (across >= 0 && !queue.contains(across) &&
                            r.faces[static_cast<std::size_t>(across)] == -1) {
                            const auto corners = adjacency_.edgeCorners(ea);
                            const int acrossCorner = corners[0] == a ? corners[1] : (corners[1] == a ? corners[0] : -1);
                            nextFaces.push_back({across, acrossCorner});
                        }
                    } else if (known != eb) {
                        fits = false;
                        break;
                    }
                }
                a = adjacency_.next(a);
                b = adjacency_.prev(b);
            }
            if (!fits) {
                continue;
            }
            r.faces[static_cast<std::size_t>(f)] = g;
            r.faces[static_cast<std::size_t>(g)] = f;
            r.side.set(static_cast<std::size_t>(f));
            r.paired += g == f ? 1 : 2;
            for (const auto& [p, q] : edgePairs) {
                if (p >= 0 && q >= 0) {
                    pair(r.edges, p, q);
                }
            }
            for (const auto& [p, q] : pointPairs) {
                pair(r.points, p, q);
            }
            for (const auto& [h, c] : nextFaces) {
                if (!queue.contains(h)) {
                    entry[static_cast<std::size_t>(h)] = c;
                    queue.pushBack(h);
                }
            }
        }
        return r;
    }

    // The winning pairs become the cache: sides, normal and point.
    SymmetryCache take(Pairing r, int axis) const
    {
        const auto pointCount = points_.size();
        SymmetryCache::Maps maps{std::move(r.points), BitArray(pointCount), std::move(r.faces), std::move(r.side)};
        const Topology& topology = mesh_.topology();
        for (int f = 0; f < faceCount_; ++f) {
            if (!maps.faceSide.test(static_cast<std::size_t>(f))) {
                continue;
            }
            for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
                maps.pointSide.setBounded(static_cast<std::size_t>(topology.vertex(c)));
            }
        }
        Vec3f d;
        Vec3f sum;
        int n = 0;
        for (std::size_t v = 0; v < pointCount; ++v) {
            const int m = maps.points[v];
            if (m < 0) {
                continue;
            }
            if (m != static_cast<int>(v) && maps.pointSide.test(v)) {
                d = d + (points_[v] - at(m));
            }
            sum = sum + points_[v];
            ++n;
        }
        normalizeInPlace(d);
        Vec3f normal;
        if (d.x > snapLimit || d.x < -snapLimit) {
            normal = {1.0f, 0.0f, 0.0f};
        } else if (d.y > snapLimit || d.y < -snapLimit) {
            normal = {0.0f, 1.0f, 0.0f};
        } else if (d.z > snapLimit || d.z < -snapLimit) {
            normal = {0.0f, 0.0f, 1.0f};
        } else if (axis < 0 || !(axisValue(d, axis) < 0.0f)) {
            normal = d;
        } else {
            normal = {-d.x, -d.y, -d.z};
        }
        if (0.0f > dot(d, normal)) {
            maps.faceSide.flipAll();
            maps.pointSide.flipAll();
        }
        // A locked plane keeps normal and point; without paired points the
        // point stays where it was.
        Vec3f planeNormal = normal;
        Vec3f planePoint;
        if (lock_) {
            planeNormal = lock_->normal;
            planePoint = lock_->point;
        } else if (n > 0) {
            const auto count = static_cast<float>(n);
            planePoint = {sum.x / count, sum.y / count, sum.z / count};
        }
        if (r.paired != faceCount_) {
            pairByDistance(maps.points, planeNormal, planePoint);
        }
        return SymmetryCache(std::move(maps), std::move(r.edges), r.paired, lock_);
    }

    // Points without partner, sorted by their distance along the normal, pair
    // with the nearest unpaired point near their reflection (closer than half
    // the shortest edge at such points; one closer than 1/400 of it is taken
    // at once). Points just above the plane pair with themselves unless a
    // partner is found.
    void pairByDistance(std::vector<int>& mirror, Vec3f normal, Vec3f planePoint) const
    {
        const float d0 = dot(planePoint, normal);
        std::vector<int> open;
        for (std::size_t v = 0; v < mirror.size(); ++v) {
            if (mirror[v] < 0) {
                open.push_back(static_cast<int>(v));
            }
        }
        float shortest = FLT_MAX;
        for (const int u : open) {
            for (const int e : adjacency_.pointEdges(u)) {
                const auto ends = adjacency_.edgeVertices(e);
                const Vec3f d = at(ends[0]) - at(ends[1]);
                const float q = dot(d, d);
                shortest = q < shortest ? q : shortest;
            }
        }
        const float h = static_cast<float>(std::sqrt(static_cast<double>(shortest)) * 0.5);
        const float strict = (h / strictPairDivisor) * (h / strictPairDivisor);
        auto key = [&](int v) { return dot(normal, at(v)); };
        std::vector<int> sorted = open;
        util::Introsort::sort(sorted.begin(), sorted.end(), [&](int a, int b) { return key(a) < key(b); });
        auto partner = [&](int v) -> int& { return mirror[static_cast<std::size_t>(v)]; };
        int j = static_cast<int>(sorted.size()) - 1;
        for (const int u : sorted) {
            if (partner(u) >= 0) {
                continue;
            }
            const float su = key(u);
            if (su > d0) {
                if (su > d0 + h) {
                    break;
                }
                partner(u) = u;
            }
            const float t = 2.0f * d0 - su;
            const float upper = t + h;
            const float lower = t - h;
            if (j < 0) {
                continue;
            }
            int first = -1;
            float best = h * h;
            int nearest = -1;
            for (int k = j; k >= 0; --k) {
                const int w = sorted[static_cast<std::size_t>(k)];
                if (partner(w) >= 0) {
                    continue;
                }
                const float sw = key(w);
                if (sw > upper) {
                    continue;
                }
                if (lower > sw) {
                    break;
                }
                if (first == -1) {
                    first = k;
                }
                const float q = (sw - d0) * 2.0f;
                const Vec3f reflected = at(w) - normal * q;
                const Vec3f diff = at(u) - reflected;
                const float dist = dot(diff, diff);
                if (strict > dist) {
                    pair(mirror, u, w);
                    break;
                }
                if (best > dist) {
                    best = dist;
                    nearest = w;
                }
            }
            if (partner(u) == -1 && nearest != -1) {
                pair(mirror, u, nearest);
            }
            if (first != -1) {
                j = first;
            }
        }
    }

    const Mesh& mesh_;
    const Adjacency& adjacency_;
    const std::vector<Vec3f>& points_;
    int faceCount_;
    float narrowLimit_;
    std::optional<SplitNormals> splitNormals_;
    BitArray onLine_;  // edges on a followed line, for the whole run
    BitArray crossed_; // edges crossed by a line, for the whole run
    std::optional<SymmetryCache::LockedPlane> lock_;
};

} // namespace

std::optional<SymmetryCache> SymmetryDetector::run(const Mesh& mesh, const SymmetryCache* previous)
{
    return Matcher(mesh, previous).run();
}

} // namespace areko::remesher::prep
