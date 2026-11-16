#include "symmetry_cut.h"

#include "edge_splitter.h"
#include "mesh_repair.h"

#include "api/status.h"
#include "mesh/modification.h"
#include "mesh/sharp_edges.h"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr double toleranceFactor = 1e-6;
constexpr double planeSnapFactor = 1e-7;
constexpr double crossSnapFactor = 0.0005;

// Plane position with the bounding box middle snapped to 0 below `limit`.
float planeOf(const SymmetryCut::Job& job, const Bounds& bounds, float limit)
{
    if (!job.local) {
        return job.usePresetPlane ? job.planePos : 0.0f;
    }
    const float c = axisValue(center(bounds), job.axis);
    return std::fabs(c) < limit ? 0.0f : c;
}

// Gives every edge from one side of the plane to the other (neither end
// within the tolerance) a point on the plane. Returns the plane position.
float splitAtPlane(Mesh& mesh, const SymmetryCut::Job& job)
{
    if (mesh.topology().mode() != FaceMode::Triangles) {
        throw api::StatusError(ARM_ERROR_INTERNAL, "the plane split needs a triangle mesh");
    }
    const int a = job.axis;
    const float tol = job.tolerance;
    const float c = planeOf(job, mesh.bounds(), tol);
    const std::vector<Vec3f>& points = mesh.points();
    // Side per point: 0 below, 1 on the plane, 2 above (also NaN).
    std::vector<unsigned char> side(points.size(), 1);
    for (std::size_t v = 0; v < points.size(); ++v) {
        const float x = axisValue(points[v], a);
        if (!(x <= c + tol)) {
            side[v] = 2;
        } else if (x < c - tol) {
            side[v] = 0;
        }
    }
    const Adjacency& adjacency = mesh.adjacency();
    std::vector<int> edges;
    std::vector<Vec3f> cuts;
    std::vector<float> ts;
    for (int e = 0; e < adjacency.edgeCount(); ++e) {
        const auto ends = adjacency.edgeVertices(e);
        const unsigned char sp = side[static_cast<std::size_t>(ends[0])];
        const unsigned char sq = side[static_cast<std::size_t>(ends[1])];
        if (sp == 1 || sq == 1 || sp == sq) {
            continue;
        }
        const Vec3f p = points[static_cast<std::size_t>(ends[0])];
        const Vec3f q = points[static_cast<std::size_t>(ends[1])];
        const float t = (c - axisValue(p, a)) / (axisValue(q, a) - axisValue(p, a));
        Vec3f cut = q * t + p * (1.0f - t);
        axisValue(cut, a) = c;
        edges.push_back(e);
        cuts.push_back(cut);
        ts.push_back(t);
    }
    EdgeSplitter::splitAt(mesh, edges, cuts, ts);
    return c;
}

// Cutting and mirroring on a copy of the points; sequential, later edges see
// the points moved before them.
class Cutter {
public:
    Cutter(Mesh& mesh, const SymmetryCut::Job& job)
        : mesh_(mesh), job_(job), adjacency_(mesh.adjacency()), topology_(mesh.topology()), pos_(mesh.points()),
          axis_(job.axis), tol_(job.tolerance)
    {
        const Bounds& bounds = mesh.bounds();
        diagonal_ = diagonal(bounds);
        c_ = planeOf(job, bounds, static_cast<float>(static_cast<double>(diagonal_) * planeSnapFactor));
        keepPointOrigin_ = job.keepPointOrigin && mesh.maps().hasKind(ElementKind::Point);
        keepFaceOrigin_ = job.keepFaceOrigin && mesh.maps().hasKind(ElementKind::Face);
        if (job.keepHardEdges) {
            hard_ = SharpEdges::halfEdgeBits(mesh);
        }
    }

    std::optional<SymmetryCache> run()
    {
        const int pointCount = mesh_.pointCount();
        classifyPoints();
        moveCrossingEnds();
        buildPoints();
        buildFaces();
        apply();
        if (!job_.buildCache || job_.cutOnly) {
            return std::nullopt;
        }
        return cache(pointCount);
    }

private:
    [[nodiscard]] float x(int v) const { return axisValue(pos_[static_cast<std::size_t>(v)], axis_); }
    float& x(int v) { return axisValue(pos_[static_cast<std::size_t>(v)], axis_); }
    [[nodiscard]] bool kept(int v) const { return newIndex_[static_cast<std::size_t>(v)] != -1; }
    [[nodiscard]] bool onPlane(int v) const { return onPlane_.test(static_cast<std::size_t>(v)); }

    // Kept points (negative side) in index order; those within the tolerance
    // snap onto the plane.
    void classifyPoints()
    {
        const auto n = pos_.size();
        onPlane_ = BitArray(n);
        newIndex_.assign(n, -1);
        for (std::size_t v = 0; v < n; ++v) {
            const float xv = axisValue(pos_[v], axis_);
            if (!(xv <= c_ + tol_)) {
                continue;
            }
            if (xv >= c_ - tol_) {
                onPlane_.set(v);
                axisValue(pos_[v], axis_) = c_;
            }
            newIndex_[v] = kept_++;
        }
    }

    // Of every half edge with exactly one kept end and no end on the plane,
    // one end moves onto the plane: the one on the positive side unless the
    // plane lies closer than the cross ratio to the negative end.
    void moveCrossingEnds()
    {
        float snap = tol_ + tol_;
        if (static_cast<double>(snap) < static_cast<double>(diagonal_) * crossSnapFactor) {
            snap = static_cast<float>(static_cast<double>(diagonal_) * crossSnapFactor);
        }
        const float r = job_.crossRatio;
        const float rest = 1.0f - r;
        for (int f = 0; f < topology_.faceCount(); ++f) {
            for (int k = topology_.faceBegin(f); k < topology_.faceEnd(f); ++k) {
                const int cur = adjacency_.vertex(k);
                const int nxt = adjacency_.vertex(adjacency_.next(k));
                if (kept(cur) == kept(nxt) || onPlane(cur) || onPlane(nxt)) {
                    continue;
                }
                const float s = (x(nxt) - c_) / (x(nxt) - x(cur));
                const int v = !(x(nxt) > x(cur)) ? (s < r ? nxt : cur) : (s > rest ? cur : nxt);
                if (std::fabs(x(v) - c_) < snap) {
                    x(v) = c_;
                } else {
                    moveAlongNeighbour(v);
                }
                onPlane_.set(static_cast<std::size_t>(v));
                if (!kept(v)) {
                    newIndex_[static_cast<std::size_t>(v)] = kept_++;
                }
            }
        }
    }

    // Distance of w to the plane when it lies on the other side of v (or on the plane).
    [[nodiscard]] bool across(int v, int w, float& dist) const
    {
        if (!(x(v) < c_)) {
            if (!(x(w) <= c_)) {
                return false;
            }
            dist = c_ - x(w);
        } else {
            if (!(x(w) >= c_)) {
                return false;
            }
            dist = x(w) - c_;
        }
        return true;
    }

    // Moves v along the edge to its most distant neighbour across the plane
    // (next point of a corner, then its previous point; when the next point is
    // on the same side, the previous one of that corner is not looked at). A
    // corner whose half edge is sharp decides at once when its next point lies
    // across. Without such a neighbour v only moves onto the plane.
    void moveAlongNeighbour(int v)
    {
        float best = -1.0f;
        int bestCorner = -1;
        for (const int cc : adjacency_.pointCorners(v)) {
            const int nextCorner = adjacency_.next(cc);
            const int prevCorner = adjacency_.prev(cc);
            float dist = 0.0f;
            const bool nextAcross = across(v, adjacency_.vertex(nextCorner), dist);
            if (hard_ && hard_->testBounded(static_cast<std::size_t>(cc)) && nextAcross) {
                bestCorner = nextCorner;
                break;
            }
            if (!nextAcross) {
                continue;
            }
            if (dist > best) {
                best = dist;
                bestCorner = nextCorner;
            }
            if (!across(v, adjacency_.vertex(prevCorner), dist)) {
                continue;
            }
            if (dist > best) {
                best = dist;
                bestCorner = prevCorner;
            }
        }
        if (bestCorner != -1) {
            const int w = adjacency_.vertex(bestCorner);
            const float lambda = (x(w) - c_) / (x(w) - x(v));
            Vec3f& p = pos_[static_cast<std::size_t>(v)];
            p = pos_[static_cast<std::size_t>(w)] * (1.0f - lambda) + p * lambda;
        }
        x(v) = c_;
    }

    // Kept points in their new order; with mirroring the copies of the points
    // off the plane follow in the old order of their originals.
    void buildPoints()
    {
        const float twoC = c_ + c_;
        mirror_.assign(static_cast<std::size_t>(kept_), 0);
        newPoints_.assign(static_cast<std::size_t>(kept_), Vec3f{});
        for (std::size_t v = 0; v < pos_.size(); ++v) {
            const int ni = newIndex_[v];
            if (ni < 0) {
                continue;
            }
            newPoints_[static_cast<std::size_t>(ni)] = pos_[v];
            if (job_.cutOnly) {
                // Keep-half: a point on the plane is its own mirror, the others have none.
                mirror_[static_cast<std::size_t>(ni)] = onPlane_.test(v) ? ni : -1;
            } else if (onPlane_.test(v)) {
                mirror_[static_cast<std::size_t>(ni)] = ni;
            } else {
                Vec3f p = pos_[v];
                axisValue(p, axis_) = twoC - axisValue(p, axis_);
                mirror_[static_cast<std::size_t>(ni)] = static_cast<int>(newPoints_.size());
                newPoints_.push_back(p);
            }
        }
    }

    // Faces with more than two kept corners, not all of them on the plane,
    // whose kept corners do not all mirror onto themselves; then, with
    // mirroring, their mirror faces in reversed corner order.
    void buildFaces()
    {
        const bool hard = hard_.has_value();
        starts_.assign(1, 0);
        for (int f = 0; f < topology_.faceCount(); ++f) {
            int keptCorners = 0;
            bool allOnPlane = true;
            for (int k = topology_.faceBegin(f); k < topology_.faceEnd(f); ++k) {
                const int v = adjacency_.vertex(k);
                keptCorners += kept(v) ? 1 : 0;
                allOnPlane = allOnPlane && onPlane(v);
            }
            if (allOnPlane || keptCorners <= 2) {
                continue;
            }
            const std::size_t begin = corners_.size();
            bool allSelf = true;
            for (int k = topology_.faceBegin(f); k < topology_.faceEnd(f); ++k) {
                const int ni = newIndex_[static_cast<std::size_t>(adjacency_.vertex(k))];
                if (ni == -1) {
                    continue;
                }
                corners_.push_back(ni);
                if (hard) {
                    cornerOrigin_.push_back(k);
                }
                allSelf = allSelf && ni == mirror_[static_cast<std::size_t>(ni)];
            }
            if (allSelf || corners_.size() - begin < 3) {
                corners_.resize(begin);
                cornerOrigin_.resize(hard ? begin : 0);
                continue;
            }
            starts_.push_back(static_cast<int>(corners_.size()));
            if (keepFaceOrigin_) {
                faceOrigin_.push_back(f);
            }
        }
        const std::size_t half = starts_.size() - 1;
        if (job_.cutOnly) {
            return;
        }
        // The half edge v_k -> v_k+1 becomes v'_k+1 -> v'_k: the mirrored corner
        // of v_k takes the origin of corner k - 1 (the first one of the last).
        for (std::size_t i = 0; i < half; ++i) {
            const auto s = static_cast<std::size_t>(starts_[i]);
            const auto e = static_cast<std::size_t>(starts_[i + 1]) - 1;
            for (std::size_t p = e + 1; p-- > s;) {
                const int m = mirror_[static_cast<std::size_t>(corners_[p])];
                corners_.push_back(m);
                if (hard) {
                    const int origin = cornerOrigin_[p == s ? e : p - 1];
                    cornerOrigin_.push_back(origin);
                }
            }
            starts_.push_back(static_cast<int>(corners_.size()));
        }
        if (keepFaceOrigin_) {
            const std::vector<int> first = faceOrigin_;
            faceOrigin_.insert(faceOrigin_.end(), first.begin(), first.end());
        }
    }

    void apply()
    {
        Modification change;
        if (keepPointOrigin_) {
            std::vector<int> origin(newPoints_.size(), -1);
            for (std::size_t v = 0; v < pos_.size(); ++v) {
                const int ni = newIndex_[v];
                if (ni == -1) {
                    continue;
                }
                origin[static_cast<std::size_t>(ni)] = static_cast<int>(v);
                const int m = mirror_[static_cast<std::size_t>(ni)];
                if (m != -1) {
                    origin[static_cast<std::size_t>(m)] = static_cast<int>(v);
                }
            }
            change.points.newToOld = std::move(origin);
        }
        // Triangles stay triangles; a quad mesh keeps its mode only while every face has four corners.
        FaceMode mode = topology_.mode();
        const int size = mode == FaceMode::Quads ? 4 : 3;
        for (std::size_t i = 0; mode != FaceMode::Polygons && i + 1 < starts_.size(); ++i) {
            if (starts_[i + 1] - starts_[i] != size) {
                mode = FaceMode::Polygons;
            }
        }
        if (mode == FaceMode::Polygons) {
            change.faceStart = std::move(starts_);
        }
        change.replacePoints = std::move(newPoints_);
        change.cornerVertex = std::move(corners_);
        change.faceMode = mode;
        if (keepFaceOrigin_) {
            change.faces.newToOld = std::move(faceOrigin_);
        }
        if (hard_) {
            change.halfEdges.newToOld = std::move(cornerOrigin_);
        }
        MeshModifier::apply(mesh_, change);
    }

    // Mirror pairs: kept points side 0, copies side 1; points on the plane
    // whose faces form more than one fan are split, and each copy pairs with
    // its source.
    SymmetryCache cache(int pointCountBefore)
    {
        const auto n1 = static_cast<std::size_t>(mesh_.pointCount());
        std::vector<int> pointMirror = mirror_;
        pointMirror.resize(n1, -1);
        BitArray pointSide(n1);
        for (std::size_t i = static_cast<std::size_t>(kept_); i < n1; ++i) {
            pointSide.set(i);
        }
        BitArray plane(static_cast<std::size_t>(pointCountBefore));
        for (std::size_t i = 0; i < static_cast<std::size_t>(kept_); ++i) {
            const int m = pointMirror[i];
            if (m == -1) {
                continue;
            }
            pointMirror[static_cast<std::size_t>(m)] = static_cast<int>(i);
            if (static_cast<std::size_t>(m) == i) {
                plane.set(i);
            }
        }
        const std::vector<int> sources = MeshRepair::splitNonManifoldPoints(mesh_, &plane);
        const auto count = static_cast<std::size_t>(mesh_.pointCount());
        if (sources.empty()) {
            pointMirror.resize(count, -1);
        } else {
            const auto before = pointMirror.size();
            pointMirror.resize(count, 0);
            for (std::size_t j = 0; j < sources.size(); ++j) {
                pointMirror[static_cast<std::size_t>(sources[j])] = static_cast<int>(before + j);
                pointMirror[before + j] = sources[j];
            }
        }
        pointSide.resize(count, false);

        const auto faces = static_cast<std::size_t>(mesh_.faceCount());
        const std::size_t h = faces / 2;
        std::vector<int> faceMirror(faces, 0);
        BitArray faceSide(faces);
        for (std::size_t i = 0; i < h; ++i) {
            faceMirror[i] = static_cast<int>(h + i);
            faceMirror[h + i] = static_cast<int>(i);
        }
        for (std::size_t i = h; i < faces; ++i) {
            faceSide.set(i);
        }
        // The plane stays locked to the cut axis; its point is not kept.
        SymmetryCache::LockedPlane locked;
        locked.axis = axis_;
        axisValue(locked.normal, axis_) = 1.0f;
        return SymmetryCache({std::move(pointMirror), std::move(pointSide), std::move(faceMirror), std::move(faceSide)},
                             locked);
    }

    Mesh& mesh_;
    const SymmetryCut::Job& job_;
    const Adjacency& adjacency_;
    const Topology& topology_;
    std::vector<Vec3f> pos_;
    int axis_;
    float tol_;
    float diagonal_ = 0.0f;
    float c_ = 0.0f;
    bool keepPointOrigin_ = false;
    bool keepFaceOrigin_ = false;
    std::optional<BitArray> hard_;
    BitArray onPlane_;
    std::vector<int> newIndex_;
    int kept_ = 0;
    std::vector<int> mirror_;
    std::vector<Vec3f> newPoints_;
    std::vector<int> starts_;
    std::vector<int> corners_;
    std::vector<int> cornerOrigin_;
    std::vector<int> faceOrigin_;
};

} // namespace

SymmetryCut::Result SymmetryCut::run(Mesh& mesh, Job job)
{
    Result result;
    if (job.tolerance < 0.0f) {
        job.tolerance = static_cast<float>(static_cast<double>(diagonal(mesh.bounds())) * toleranceFactor);
    }
    if (job.splitEdgesFirst) {
        result.planePos = splitAtPlane(mesh, job);
    }
    result.cache = Cutter(mesh, job).run();
    return result;
}

} // namespace areko::remesher::prep
