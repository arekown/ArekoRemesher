#include "flat_triangle_check.h"

#include "mesh/measures.h"
#include "mesh/modification.h"
#include "util/introsort.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr double oneDegree = std::bit_cast<double>(0x3f91df46a2529d39ull);
constexpr double twoDegrees = std::bit_cast<double>(0x3fa1df46a2529d39ull);
constexpr float tinyAngle = std::bit_cast<float>(0x3c0efa35u);      // 0.5 degrees
constexpr float requiredFold = std::bit_cast<float>(0x3fdf66f3u);   // 100 degrees
constexpr int stageCount = 3;
constexpr int conflictRounds = 5;

// (a < b) ? a : b, so a NaN in `a` drops out and a NaN in `b` stays.
float minOf(float a, float b) { return a < b ? a : b; }

// Turns the diagonal of the two triangles at every listed edge (triangle
// meshes). Conflicts are resolved in list order first: an edge whose new
// diagonal already exists waits for that edge, an edge whose flip would
// duplicate another edge switches that one off.
class EdgeFlipper {
public:
    static void flip(Mesh& mesh, const std::vector<int>& edges)
    {
        if (edges.empty()) {
            return;
        }
        const Adjacency& adjacency = mesh.adjacency();
        const auto edgeCount = static_cast<std::size_t>(adjacency.edgeCount());
        BitArray active(edgeCount);
        for (const int e : edges) {
            active.set(static_cast<std::size_t>(e));
        }
        resolveConflicts(adjacency, edges, active);

        const std::vector<int>& old = mesh.topology().cornerVertex();
        std::vector<int> corners = old;
        std::vector<int> origin(old.size());
        std::iota(origin.begin(), origin.end(), 0);
        for (const int e : edges) {
            const auto ends = adjacency.edgeCorners(e);
            if (!active.test(static_cast<std::size_t>(e)) || ends.size() != 2) {
                continue;
            }
            const int k0 = ends[0];
            const int k1 = ends[1];
            corners[static_cast<std::size_t>(k0)] = old[static_cast<std::size_t>(adjacency.prev(k1))];
            corners[static_cast<std::size_t>(k1)] = old[static_cast<std::size_t>(adjacency.prev(k0))];
            const int faces[2] = {k0 / 3, k1 / 3};
            for (const int f : faces) {
                for (int c = 3 * f; c < 3 * f + 3; ++c) {
                    origin[static_cast<std::size_t>(c)] =
                        oldHalfEdge(adjacency, old, faces, corners[static_cast<std::size_t>(c)],
                                    corners[static_cast<std::size_t>(adjacency.next(c))]);
                }
            }
        }

        Modification change;
        change.cornerVertex = std::move(corners);
        change.faceMode = FaceMode::Triangles;
        change.points.state = 1;
        change.faces.state = 1;
        change.halfEdges.newToOld = std::move(origin);
        MeshModifier::apply(mesh, change);
    }

private:
    // First corner of the two old faces whose half edge runs from u to v, or -1.
    static int oldHalfEdge(const Adjacency& adjacency, const std::vector<int>& old, const int (&faces)[2], int u, int v)
    {
        for (const int f : faces) {
            for (int c = 3 * f; c < 3 * f + 3; ++c) {
                if (old[static_cast<std::size_t>(c)] == u &&
                    old[static_cast<std::size_t>(adjacency.next(c))] == v) {
                    return c;
                }
            }
        }
        return -1;
    }

    static void resolveConflicts(const Adjacency& adjacency, const std::vector<int>& edges, BitArray& active)
    {
        const auto edgeCount = static_cast<std::size_t>(adjacency.edgeCount());
        BitArray finished(edgeCount);
        std::vector<std::uint32_t> stamp(edgeCount, 0);
        std::uint32_t counter = 1;
        std::vector<int> list = edges;
        for (int round = 0;;) {
            std::vector<int> deferred;
            for (const int e : list) {
                const auto ends = adjacency.edgeCorners(e);
                if (ends.size() != 2) {
                    continue;
                }
                const int a = adjacency.vertex(adjacency.prev(ends[0]));
                const int b = adjacency.vertex(adjacency.prev(ends[1]));
                const int diagonal = adjacency.edgeBetween(a, b);
                if (diagonal == -1) {
                    // Edges opposite to both a and b would appear twice after the flip.
                    if (++counter == 0) {
                        std::fill(stamp.begin(), stamp.end(), 0u);
                        counter = 1;
                    }
                    for (const int q : adjacency.pointCorners(a)) {
                        const int x = adjacency.cornerEdge(adjacency.next(q));
                        if (x != -1) {
                            stamp[static_cast<std::size_t>(x)] = counter;
                        }
                    }
                    for (const int q : adjacency.pointCorners(b)) {
                        const int x = adjacency.cornerEdge(adjacency.next(q));
                        if (x != -1 && x != e && stamp[static_cast<std::size_t>(x)] == counter) {
                            active.set(static_cast<std::size_t>(x), false);
                        }
                    }
                    finished.set(static_cast<std::size_t>(e));
                } else if (active.test(static_cast<std::size_t>(diagonal))) {
                    if (finished.test(static_cast<std::size_t>(diagonal))) {
                        finished.set(static_cast<std::size_t>(e));
                    } else {
                        deferred.push_back(e);
                    }
                } else {
                    finished.set(static_cast<std::size_t>(e));
                    active.set(static_cast<std::size_t>(e), false);
                }
            }
            if (deferred.empty()) {
                // A round after the first that defers nothing switches its whole list off.
                if (round > 0) {
                    switchOff(list, active);
                }
                return;
            }
            ++round;
            list = std::move(deferred);
            if (round == conflictRounds) {
                switchOff(list, active);
                return;
            }
        }
    }

    static void switchOff(const std::vector<int>& list, BitArray& active)
    {
        for (const int e : list) {
            active.set(static_cast<std::size_t>(e), false);
        }
    }
};

// Chooses the edges to flip at acute corners of one triangle mesh; in mode 1
// required edges come first.
class FlipChoice {
public:
    FlipChoice(const Mesh& mesh, const BitArray& locked, const FlatTriangleBlock& block)
        : adjacency_(mesh.adjacency()), points_(mesh.points()), angles_(cornerAngles(mesh)),
          dihedral_(mesh.dihedralAngles()), normals_(mesh.triangleGeometry().normals()), locked_(locked),
          required_(requiredEdges(block)), angleMain_(block.angleMain), angleSmall_(block.angleSmall)
    {
        const TriangleGeometry& geometry = mesh.triangleGeometry();
        const float limit = geometry.smallAreaLimit();
        const auto& doubleAreas = geometry.doubleAreas();
        flat_ = BitArray(static_cast<std::size_t>(mesh.faceCount()));
        for (int f = 0; f < mesh.faceCount(); ++f) {
            if (doubleAreas[static_cast<std::size_t>(f)] > limit) {
                continue;
            }
            for (int c = 3 * f; c < 3 * f + 3; ++c) {
                if (angles_[static_cast<std::size_t>(c)] < tinyAngle) {
                    flat_.set(static_cast<std::size_t>(f));
                }
            }
        }
    }

    // Smallest corner angle up to the main angle before any flip (FLT_MAX
    // without such a corner, -1 without corners).
    [[nodiscard]] float minAngle() const noexcept { return minAngle_; }

    [[nodiscard]] std::vector<int> choose()
    {
        std::vector<int> chosen;
        // Required edges are flipped without any check unless they are locked
        // already; their neighbour edges become locked.
        for (const int e : required_) {
            if (locked_.testBounded(static_cast<std::size_t>(e))) {
                continue;
            }
            lock(e);
            chosen.push_back(e);
            lockNeighbours(e);
        }

        const float s1 = angleMain_ / 3.0f;
        const float s2 = (angleMain_ + angleMain_) / 3.0f;
        const float s3 = (angleMain_ * 3.0f) / 3.0f;
        struct Stage {
            std::vector<int> corners;
            std::vector<float> angles;
        };
        Stage stages[stageCount];
        const auto cornerCount = static_cast<int>(angles_.size());
        minAngle_ = cornerCount == 0 ? -1.0f : std::numeric_limits<float>::max();
        for (int c = 0; c < cornerCount; ++c) {
            const float a = angles_[static_cast<std::size_t>(c)];
            if (a > s3) {
                continue;
            }
            Stage& stage = stages[a > s2 ? 2 : (a > s1 ? 1 : 0)];
            stage.corners.push_back(c);
            stage.angles.push_back(a);
            if (a < minAngle_) {
                minAngle_ = a;
            }
        }

        for (const Stage& stage : stages) {
            std::vector<int> order(stage.corners.size());
            std::iota(order.begin(), order.end(), 0);
            util::Introsort::sort(order.begin(), order.end(), [&stage](int a, int b) {
                return stage.angles[static_cast<std::size_t>(a)] < stage.angles[static_cast<std::size_t>(b)];
            });
            for (const int i : order) {
                const int choice = chooseAt(stage.corners[static_cast<std::size_t>(i)],
                                            stage.angles[static_cast<std::size_t>(i)]);
                if (choice >= 0) {
                    chosen.push_back(choice);
                }
            }
        }
        return chosen;
    }

private:
    [[nodiscard]] bool isLocked(int e) const { return locked_.testSigned(e); }
    void lock(int e) { locked_.setSigned(e); }

    // Locks the other edges of the triangles at e.
    void lockNeighbours(int e)
    {
        for (const int k : adjacency_.edgeCorners(e)) {
            lock(adjacency_.cornerEdge(adjacency_.next(k)));
            lock(adjacency_.cornerEdge(adjacency_.prev(k)));
        }
    }

    // Required edges of mode 1 with a point range (begin >= 0), ascending:
    // folded by at least 100 degrees (or NaN), both end points in the range.
    [[nodiscard]] std::vector<int> requiredEdges(const FlatTriangleBlock& block) const
    {
        std::vector<int> out;
        if (block.mode != 1 || block.requiredBegin < 0) {
            return out;
        }
        auto inRange = [&block](int v) { return block.requiredBegin <= v && v < block.requiredEnd; };
        for (int e = 0; e < adjacency_.edgeCount(); ++e) {
            if (requiredFold > dihedral_[static_cast<std::size_t>(e)]) {
                continue;
            }
            const auto ends = adjacency_.edgeVertices(e);
            if (inRange(ends[0]) && inRange(ends[1])) {
                out.push_back(e);
            }
        }
        return out;
    }

    // Edge to flip for one acute corner, or -1. The chosen edge and the four
    // other edges of its triangles become locked at once.
    int chooseAt(int c, float threshold)
    {
        const int e1 = adjacency_.cornerEdge(c);
        const int e2 = adjacency_.cornerEdge(adjacency_.prev(c));
        if (e1 == -1 || e2 == -1) {
            if (e1 != -1) {
                lock(e1);
            }
            if (e2 != -1) {
                lock(e2);
            }
            return -1;
        }
        const bool flatFace = flat_.test(static_cast<std::size_t>(c / 3));
        float q1 = 0.0f;
        float q2 = 0.0f;
        const bool ok1 = !isLocked(e1) && (dihedral_[static_cast<std::size_t>(e1)] <= angleSmall_ || flatFace) &&
                         evaluate(threshold, e1, q1);
        const bool ok2 = !isLocked(e2) && (dihedral_[static_cast<std::size_t>(e2)] <= angleSmall_ || flatFace) &&
                         evaluate(threshold, e2, q2);
        int choice = -1;
        if (ok1 && ok2) {
            choice = q2 > q1 ? e2 : e1;
        } else if (ok1) {
            choice = e1;
        } else if (ok2) {
            choice = e2;
        } else {
            return -1;
        }
        lockNeighbours(choice);
        lock(choice);
        return choice;
    }

    [[nodiscard]] Vec3f point(int corner) const
    {
        return points_[static_cast<std::size_t>(adjacency_.vertex(corner))];
    }

    // Accepts the flip of e when the quad stays convex towards the triangle
    // normals, no new corner at the far points is smaller than `threshold`
    // and the new fold stays within twice the old one (2 degrees for an
    // almost flat edge). `quality` is the smallest of the four new angles.
    bool evaluate(float threshold, int e, float& quality) const
    {
        const auto ends = adjacency_.edgeCorners(e);
        if (ends.size() != 2) {
            return false;
        }
        const int k0 = ends[0];
        const int k1 = ends[1];
        const Vec3f p1 = point(k0);
        const Vec3f p2 = point(k1);
        const Vec3f p3 = point(adjacency_.prev(k0));
        const Vec3f p4 = point(adjacency_.prev(k1));
        const int f0 = k0 / 3;
        const int f1 = k1 / 3;
        Vec3f n0 = normals_[static_cast<std::size_t>(f0)];
        Vec3f n1 = normals_[static_cast<std::size_t>(f1)];
        const bool flat0 = flat_.test(static_cast<std::size_t>(f0));
        const bool flat1 = flat_.test(static_cast<std::size_t>(f1));
        if (flat0 && flat1) {
            return false;
        }
        if (flat0) {
            n0 = n1;
        }
        if (flat1) {
            n1 = n0;
        }
        const Vec3f d = p4 - p3;
        const Vec3f u1 = p1 - p3;
        const Vec3f u2 = p2 - p3;
        const float a1 = signedVectorAngle(u1, d, cross(u1, d), n0);
        const float a2 = signedVectorAngle(u2, d, cross(d, u2), n0);
        if (a1 >= piF || a2 >= piF) {
            return false;
        }
        const Vec3f u3 = p4 - p2;
        const Vec3f u4 = p4 - p1;
        const float a3 = signedVectorAngle(u3, d, cross(u3, d), n1);
        const float a4 = signedVectorAngle(d, u4, cross(d, u4), n1);
        if (a3 >= piF || a4 >= piF) {
            return false;
        }
        const float smallest = minOf(a4, minOf(a3, minOf(a2, a1)));
        if (threshold > smallest) {
            return false;
        }
        Vec3f m1 = cross(p4 - p1, p3 - p1);
        normalizeInPlace(m1);
        Vec3f m2 = cross(p3 - p2, p4 - p2);
        normalizeInPlace(m2);
        const float fold = static_cast<float>(angleOfCosine(dot(m1, m2)));
        const float before = dihedral_[static_cast<std::size_t>(e)];
        const bool ok = static_cast<double>(before) < oneDegree ? !(static_cast<double>(fold) > twoDegrees)
                                                                 : !(fold > before + before);
        if (ok) {
            quality = smallest;
        }
        return ok;
    }

    const Adjacency& adjacency_;
    const std::vector<Vec3f>& points_;
    std::vector<float> angles_;
    const std::vector<float>& dihedral_;
    const std::vector<Vec3f>& normals_;
    BitArray locked_;
    std::vector<int> required_;
    BitArray flat_;
    float angleMain_;
    float angleSmall_;
    float minAngle_ = -1.0f;
};

} // namespace

void FlatTriangleCheck::run(Mesh& mesh, const BitArray& locked, FlatTriangleBlock& block)
{
    block.count = 0;
    block.again = false;
    if (mesh.topology().mode() != FaceMode::Triangles) {
        return;
    }
    std::vector<int> chosen;
    float minAngle = -1.0f;
    {
        FlipChoice choice(mesh, locked, block);
        chosen = choice.choose();
        minAngle = choice.minAngle();
    }
    block.count = static_cast<int>(chosen.size());
    EdgeFlipper::flip(mesh, chosen);
    // After a flip another run follows when a corner below 1 degree was left
    // before it, or by round: mode 0 up to round 3, mode 3 in round 0; modes
    // 1 and 2 never.
    if (block.count != 0) {
        const bool tiny = static_cast<double>(minAngle) < oneDegree;
        if (block.mode == 0) {
            block.again = tiny || block.round <= 3;
        } else if (block.mode == 3) {
            block.again = tiny || block.round <= 0;
        }
    }
}

} // namespace areko::remesher::prep
