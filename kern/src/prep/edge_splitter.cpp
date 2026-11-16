#include "edge_splitter.h"

#include "symmetry_detector.h"

#include "mesh/measures.h"
#include "mesh/modification.h"

#include <algorithm>
#include <array>
#include <bit>
#include <climits>
#include <cstddef>
#include <initializer_list>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr float sharpCornerLimit = std::bit_cast<float>(0x3e567750u);  // 12 degrees
constexpr float foldLimit = std::bit_cast<float>(0x3db2b8c2u);         // 5 degrees
constexpr double piOfFloat = std::bit_cast<double>(0x400921fb60000000ull);
constexpr double angleMargin = std::bit_cast<double>(0x3eb0c6f7a0b5ed8dull); // 1e-6

// New points of one pass and the rebuild of the triangles around them.
class SplitPass {
public:
    // `mirror` (may be null): mirror pairs of the mesh with edge pairs.
    SplitPass(const Mesh& mesh, int rest, SymmetryCache* mirror)
        : mesh_(mesh), adjacency_(mesh.adjacency()), points_(mesh.points()), firstNew_(mesh.pointCount()),
          rest_(rest), newVertex_(static_cast<std::size_t>(adjacency_.edgeCount()), -1)
    {
        if (mirror != nullptr && !mirror->edgeMirror(adjacency_).empty()) {
            mirror_ = mirror;
        }
    }

    [[nodiscard]] int added() const noexcept { return static_cast<int>(newPoints_.size()); }
    [[nodiscard]] int end() const noexcept { return firstNew_ + added(); }
    [[nodiscard]] bool budgetReached() const noexcept { return rest_ <= added(); }

    // Phase A: edges longer than the target length and their mirror edges
    // (those without the border protection); a marked edge is not checked again.
    [[nodiscard]] BitArray markLongEdges(float targetLength, bool protectBorders) const
    {
        const float limit = targetLength * targetLength;
        BitArray marked(static_cast<std::size_t>(adjacency_.edgeCount()));
        const std::vector<int>* edgeMirror = mirror_ != nullptr ? &mirror_->edgeMirror(adjacency_) : nullptr;
        for (int e = 0; e < adjacency_.edgeCount(); ++e) {
            if (marked.test(static_cast<std::size_t>(e))) {
                continue;
            }
            if (protectBorders && adjacency_.edgeCorners(e).size() == 1) {
                continue;
            }
            const auto ends = adjacency_.edgeVertices(e);
            const Vec3f d = at(ends[0]) - at(ends[1]);
            if (limit < dot(d, d)) {
                marked.set(static_cast<std::size_t>(e));
                if (edgeMirror != nullptr && (*edgeMirror)[static_cast<std::size_t>(e)] != -1) {
                    marked.set(static_cast<std::size_t>((*edgeMirror)[static_cast<std::size_t>(e)]));
                }
            }
        }
        return marked;
    }

    // Phase B: one point for both marked edges of a sharp corner. Stops when
    // the budget is reached.
    void splitSharpCorners(const BitArray& marked, const BitArray* locked)
    {
        const std::vector<float> angles = cornerAngles(mesh_);
        const std::vector<float>& dihedral = mesh_.dihedralAngles();
        auto folded = [&](int e) {
            return dihedral[static_cast<std::size_t>(e)] >= foldLimit || (locked != nullptr && locked->testSigned(e));
        };
        for (int c = 0; c < static_cast<int>(angles.size()); ++c) {
            const float angle = angles[static_cast<std::size_t>(c)];
            if (angle > sharpCornerLimit) {
                continue;
            }
            const int pc = adjacency_.prev(c);
            const int nc = adjacency_.next(c);
            const int eN = adjacency_.cornerEdge(c);
            const int eP = adjacency_.cornerEdge(pc);
            const int eO = adjacency_.cornerEdge(nc);
            if (!(marked.testSigned(eN) && marked.testSigned(eP) && !marked.testSigned(eO)) || eN < 0 || eP < 0) {
                continue;
            }
            if (newVertex(eN) != -1 || newVertex(eP) != -1) {
                continue;
            }
            const int vc = adjacency_.vertex(c);
            const int vn = adjacency_.vertex(nc);
            const int vp = adjacency_.vertex(pc);
            if (valenceForbids(vc, eN, eP, marked)) {
                continue;
            }
            const bool foldP = folded(eP);
            const bool foldN = folded(eN);
            const auto cornersP = static_cast<int>(adjacency_.edgeCorners(eP).size());
            const auto cornersN = static_cast<int>(adjacency_.edgeCorners(eN).size());
            if (cornersP == 1 || foldP) {
                continue;
            }
            // On a border edge eN or a folded one the point lies in the middle of eN.
            const bool onEdge = cornersN == 1 || foldN;
            const Vec3f m = onEdge ? (at(vc) + at(vn)) * 0.5f : (((at(vc) + at(vc)) + at(vn)) + at(vp)) * 0.25f;
            const auto limit = static_cast<float>(piOfFloat - (static_cast<double>(angle) + angleMargin));
            if (cornersN == 2 && !anglesFit(c, eN, m, at(vc), at(vn), limit)) {
                continue;
            }
            if (cornersP == 2 && !anglesFit(pc, eP, m, at(vp), at(vc), limit)) {
                continue;
            }
            const int i = end();
            newVertex_[static_cast<std::size_t>(eN)] = i;
            newVertex_[static_cast<std::size_t>(eP)] = i;
            // The map values follow the position: (vc + vn) / 2 or (2 vc + vn + vp) / 4.
            if (onEdge) {
                addPoint(m, {vc, vn});
            } else {
                addPoint(m, {vc, vc, vn, vp});
            }
            if (budgetReached()) {
                return;
            }
        }
    }

    // Phase C: midpoints of the remaining marked edges, ascending.
    void splitMidpoints(const BitArray& marked)
    {
        for (std::ptrdiff_t e = marked.findNext(0); e >= 0; e = marked.findNext(static_cast<std::size_t>(e) + 1)) {
            if (newVertex(static_cast<int>(e)) != -1) {
                continue;
            }
            newVertex_[static_cast<std::size_t>(e)] = end();
            const auto ends = adjacency_.edgeVertices(static_cast<int>(e));
            addPoint((at(ends[0]) + at(ends[1])) * 0.5f, {ends[0], ends[1]});
            if (budgetReached()) {
                return;
            }
        }
    }

    // Given edges at given points (parameter mode): the map values come from
    // the end points of the edge with the weights 1 - t and t.
    void splitGiven(std::span<const int> edges, std::span<const Vec3f> points, std::span<const float> t)
    {
        for (std::size_t i = 0; i < edges.size(); ++i) {
            const int e = edges[i];
            if (newVertex(e) != -1) {
                continue;
            }
            newVertex_[static_cast<std::size_t>(e)] = end();
            const auto ends = adjacency_.edgeVertices(e);
            addPoint(points[i], {ends[0], ends[1]});
            weights_.insert(weights_.end(), {1.0f - t[i], t[i]});
        }
    }

    // Given edges at their midpoints, each followed by its mirror edge (mode
    // "middle"); an edge already split keeps its point.
    void splitMidpointsOf(std::span<const int> edges)
    {
        const std::vector<int>* edgeMirror = mirror_ != nullptr ? &mirror_->edgeMirror(adjacency_) : nullptr;
        auto splitOne = [&](int e) {
            newVertex_[static_cast<std::size_t>(e)] = end();
            const auto ends = adjacency_.edgeVertices(e);
            addPoint((at(ends[0]) + at(ends[1])) * 0.5f, {ends[0], ends[1]});
        };
        for (const int e : edges) {
            if (newVertex(e) != -1) {
                continue;
            }
            splitOne(e);
            if (edgeMirror != nullptr) {
                const int m = (*edgeMirror)[static_cast<std::size_t>(e)];
                if (m != -1 && newVertex(m) == -1) {
                    splitOne(m);
                }
            }
        }
    }

    // New triangle list in the order of the old faces; applies it with the new
    // points. Two split edges take the shorter diagonal or, without
    // `shortestDiagonal`, the one from the first split point (from the second
    // one when the mirror pairs put the triangle wholly on side 1). A triangle
    // without split edge that `splitFaces` marks fans out from its centroid.
    void rebuild(Mesh& mesh, bool shortestDiagonal, const BitArray* splitFaces = nullptr)
    {
        const std::vector<int>& corners = mesh.topology().cornerVertex();
        std::vector<int> out;
        std::vector<int> origin;
        out.reserve(corners.size() + 6 * newPoints_.size());
        origin.reserve(out.capacity() / 3);
        for (int f = 0; f < mesh.faceCount(); ++f) {
            auto emit = [&](int a, int b, int c) {
                out.insert(out.end(), {a, b, c});
                origin.push_back(f);
            };
            const std::array<int, 3> v{corners[3 * static_cast<std::size_t>(f)],
                                       corners[3 * static_cast<std::size_t>(f) + 1],
                                       corners[3 * static_cast<std::size_t>(f) + 2]};
            const std::array<int, 3> e{adjacency_.cornerEdge(3 * f), adjacency_.cornerEdge(3 * f + 1),
                                       adjacency_.cornerEdge(3 * f + 2)};
            if (e[0] < 0 || e[1] < 0 || e[2] < 0) {
                // A corner without edge: the triangle stays as it is.
                emit(v[0], v[1], v[2]);
                continue;
            }
            const std::array<int, 3> m{newVertex(e[0]), newVertex(e[1]), newVertex(e[2])};
            const int split = static_cast<int>(std::count_if(m.begin(), m.end(), [](int x) { return x > 0; }));
            if (split == 0) {
                if (splitFaces != nullptr && splitFaces->testBounded(static_cast<std::size_t>(f))) {
                    const Vec3f sum = (at(v[0]) + at(v[1])) + at(v[2]);
                    const int centre = end();
                    addPoint({sum.x / 3.0f, sum.y / 3.0f, sum.z / 3.0f}, {v[0], v[1], v[2]});
                    if (!weights_.empty()) {
                        weights_.insert(weights_.end(), {0.33333334f, 0.33333334f, 0.33333334f});
                    }
                    emit(v[0], v[1], centre);
                    emit(v[1], v[2], centre);
                    emit(v[2], v[0], centre);
                } else {
                    emit(v[0], v[1], v[2]);
                }
            } else if (split == 1) {
                const std::size_t k = m[0] > 0 ? 0 : (m[1] > 0 ? 1 : 2);
                emit(m[k], v[(k + 1) % 3], v[(k + 2) % 3]);
                emit(m[k], v[(k + 2) % 3], v[k]);
            } else if (split == 2) {
                // u: the unsplit edge; i the apex between the two split edges.
                const std::size_t u = m[0] > 0 ? (m[1] > 0 ? 2 : 1) : 0;
                const std::size_t i = (u + 2) % 3;
                const std::size_t l = (u + 1) % 3;
                const int mi = m[i];
                const int ml = m[l];
                const int vi = v[i];
                const int vj = v[u];
                const int vl = v[l];
                if (mi != ml) {
                    emit(vi, mi, ml);
                }
                if (mi == ml) {
                    emit(mi, vj, vl);
                } else if (!shortestDiagonal && mirror_ != nullptr && mirror_->sideOfFace(mesh.topology(), f) == 1) {
                    emit(mi, vj, ml);
                    emit(ml, vj, vl);
                } else if (!shortestDiagonal || firstDiagonalShorter(mi, ml, vj, vl)) {
                    emit(mi, vj, vl);
                    emit(mi, vl, ml);
                } else {
                    emit(ml, mi, vj);
                    emit(ml, vj, vl);
                }
            } else {
                if (m[0] != m[2]) {
                    emit(v[0], m[0], m[2]);
                }
                if (m[0] != m[1]) {
                    emit(v[1], m[1], m[0]);
                }
                if (m[2] != m[1]) {
                    emit(v[2], m[2], m[1]);
                    if (m[0] != m[1] && m[0] != m[2]) {
                        emit(m[0], m[1], m[2]);
                    }
                }
            }
        }

        Modification change;
        change.appendPoints = std::move(newPoints_);
        change.pointSources = std::move(sources_);
        if (!weights_.empty()) {
            change.pointWeights = std::move(weights_);
        }
        change.cornerVertex = std::move(out);
        change.faceMode = FaceMode::Triangles;
        change.points.state = 2; // old points keep their numbers
        change.faces.newToOld = std::move(origin);
        MeshModifier::apply(mesh, change);
    }

private:
    [[nodiscard]] Vec3f at(int v) const { return points_[static_cast<std::size_t>(v)]; }
    [[nodiscard]] Vec3f newPoint(int i) const { return newPoints_[static_cast<std::size_t>(i - firstNew_)]; }
    [[nodiscard]] int newVertex(int e) const { return newVertex_[static_cast<std::size_t>(e)]; }

    // Whether the diagonal mi - vl is not longer than ml - vj.
    [[nodiscard]] bool firstDiagonalShorter(int mi, int ml, int vj, int vl) const
    {
        const Vec3f di = newPoint(mi) - at(vl);
        const Vec3f dl = newPoint(ml) - at(vj);
        return dot(di, di) <= dot(dl, dl);
    }

    void addPoint(Vec3f p, std::initializer_list<int> from)
    {
        newPoints_.push_back(p);
        sources_.entries.insert(sources_.entries.end(), from);
        sources_.starts.push_back(static_cast<int>(sources_.entries.size()));
    }

    // Points with three or four corners, as many edges as corners, keep one
    // edge unsplit (three corners) or one point (four corners) on the other edges.
    [[nodiscard]] bool valenceForbids(int vc, int eN, int eP, const BitArray& marked) const
    {
        const auto edges = adjacency_.pointEdges(vc);
        const auto ne = static_cast<int>(adjacency_.pointCorners(vc).size());
        if (static_cast<int>(edges.size()) != ne || ne > 4 || ne <= 0) {
            return false;
        }
        int count = 0;
        std::vector<int> values;
        for (const int k : edges) {
            if (k != eN && k != eP && marked.test(static_cast<std::size_t>(k))) {
                ++count;
                const int value = newVertex(k);
                if (std::find(values.begin(), values.end(), value) == values.end()) {
                    values.push_back(value);
                }
            }
        }
        return (ne == 3 && count == 1) || (ne == 4 && count == 2 && values.size() == 1);
    }

    // In the triangle across edge e (seen from corner k, from `first` to
    // `second`), the new point m must not make an angle at the far point above
    // the limit; without a second corner of e at the other end no check.
    [[nodiscard]] bool anglesFit(int k, int e, Vec3f m, Vec3f first, Vec3f second, float limit) const
    {
        const auto ends = adjacency_.edgeCorners(e);
        const int k0 = ends.front();
        const int k1 = ends.back();
        const int o = k == k0 ? k1 : (k == k1 ? k0 : -1);
        if (o == -1) {
            return true;
        }
        const Vec3f q = at(adjacency_.vertex(adjacency_.prev(o)));
        const Vec3f n = mesh_.triangleGeometry().normals()[static_cast<std::size_t>(o / 3)];
        const Vec3f u = m - q;
        const Vec3f w1 = first - q;
        const Vec3f w2 = second - q;
        const float a1 = signedVectorAngle(u, w1, cross(u, w1), n);
        const float a2 = signedVectorAngle(u, w2, cross(w2, u), n);
        return !(limit < a1 || limit < a2);
    }

    const Mesh& mesh_;
    const Adjacency& adjacency_;
    const std::vector<Vec3f>& points_;
    SymmetryCache* mirror_ = nullptr;
    int firstNew_;
    int rest_;
    std::vector<int> newVertex_;
    std::vector<Vec3f> newPoints_;
    Csr sources_;
    std::vector<float> weights_; // parameter mode only
};

} // namespace

void EdgeSplitter::begin() noexcept
{
    passCount_ = 0;
    insertedTotal_ = 0;
    passInserted_ = 0;
}

bool EdgeSplitter::needsAnotherPass() const noexcept
{
    if (passCount_ < 1 || passInserted_ != 0) {
        if (maxPasses_ <= passCount_) {
            return maxPasses_ == -1;
        }
        return true;
    }
    return false;
}

void EdgeSplitter::pass(Mesh& mesh, float targetLength, const BitArray* locked)
{
    ++passCount_;
    const int pointCount = mesh.pointCount();
    const int rest = vertexBudget_ < 0 ? INT_MAX : vertexBudget_ - pointCount;
    if (vertexBudget_ >= 0 && rest < 1) {
        // No point fits any more: the pass reports nothing new.
        passInserted_ = 0;
        passFirstNew_ = 0;
        passSharpEnd_ = 0;
        return;
    }
    // The mirror pairs belong to the mesh as it is before this pass.
    std::optional<SymmetryCache> mirror = mirrorSplit_ ? SymmetryDetector::run(mesh, nullptr) : std::nullopt;
    SplitPass split(mesh, rest, mirror ? &*mirror : nullptr);
    const BitArray marked = split.markLongEdges(targetLength, protectBorders_);
    passFirstNew_ = pointCount;
    split.splitSharpCorners(marked, locked);
    passSharpEnd_ = split.end();
    // The midpoint phase adds its first point before it checks the budget.
    split.splitMidpoints(marked);
    passInserted_ = split.added();
    if (passInserted_ == 0) {
        return;
    }
    insertedTotal_ += passInserted_;
    split.rebuild(mesh, shortestDiagonal_);
}

void EdgeSplitter::splitMarked(Mesh& mesh, std::span<const int> edges, const BitArray* faces, SymmetryCache* mirror,
                               bool shortestDiagonal)
{
    SplitPass split(mesh, INT_MAX, mirror);
    split.splitMidpointsOf(edges);
    split.rebuild(mesh, shortestDiagonal, faces);
}

void EdgeSplitter::splitAt(Mesh& mesh, std::span<const int> edges, std::span<const Vec3f> points,
                           std::span<const float> t)
{
    SplitPass split(mesh, INT_MAX, nullptr);
    split.splitGiven(edges, points, t);
    split.rebuild(mesh, true);
}

} // namespace areko::remesher::prep
