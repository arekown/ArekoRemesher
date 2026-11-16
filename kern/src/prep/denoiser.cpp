#include "denoiser.h"

#include "face_borders.h"
#include "feature_chains.h"
#include "smoother.h"

#include "mesh/attribute_map.h"
#include "mesh/attribute_maps.h"
#include "mesh/bit_array.h"
#include "mesh/measures.h"
#include "mesh/sharp_edges.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <optional>
#include <span>
#include <vector>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr float flipLimit = 0.5f;

// What differs between the smoothing of the denoise step and the one after
// the crease detection.
struct ModeRules {
    int rounds;
    bool detection;       // weight rules and weight smoothing of the detection step
    float weightCap;
    bool zeroNonManifold; // weight 0 at the end points of non-manifold edges
};

ModeRules rulesFor(Denoiser::Mode mode)
{
    if (mode == Denoiser::Mode::AfterDetection) {
        return {4, true, 0.95f, false};
    }
    return {2, false, 0.6f, true};
}

// End points of face group borders (without the mesh border) and of the sharp
// edges (either form).
BitArray protectedPoints(const Mesh& mesh)
{
    BitArray points;
    const Adjacency& adjacency = mesh.adjacency();
    const int n = mesh.pointCount();
    const AttributeMaps& maps = mesh.maps();
    if (PrepState::useGroupBorders) {
        if (const auto* groups = maps.findAs<TypedMap<int>>(map_names::faceGroups, ElementKind::Face)) {
            addEdgeEndPoints(adjacency, FaceBorders::valueBorders(adjacency, groups->values(), false), n, points);
        }
    }
    if (PrepState::useEdgeAttribute) {
        if (const auto edges = SharpEdges::edgeBits(mesh)) {
            addEdgeEndPoints(adjacency, *edges, n, points);
        }
    }
    return points;
}

// Piecewise linear weight of a signed angle: 0 up to 1 degree, 0.2 at 2,
// 1 at 30, 2 at 90, 3 at 120 degrees, above that rising by 2 / pi per radian.
float ramp(float x)
{
    static constexpr float xs[5] = {degF(1), degF(2), degF(30), degF(90), degF(120)};
    static constexpr float ys[5] = {0.0f, 0.2f, 1.0f, 2.0f, 3.0f};
    const float slope = std::bit_cast<float>(0x3f22f983u);
    if (x <= xs[0]) {
        return (x - xs[0]) * 0.0f + ys[0];
    }
    int i = 1;
    while (xs[i] < x) {
        if (++i == 5) {
            return (x - xs[4]) * slope + ys[4];
        }
    }
    return ((x - xs[i]) * (ys[i - 1] - ys[i])) / (xs[i - 1] - xs[i]) + ys[i];
}

// Triangles whose doubled area is not above the small area limit.
BitArray degenerateTriangles(const TriangleGeometry& geometry)
{
    const std::vector<float>& doubleAreas = geometry.doubleAreas();
    const float limit = geometry.smallAreaLimit();
    BitArray out(doubleAreas.size());
    for (std::size_t t = 0; t < doubleAreas.size(); ++t) {
        if (!(doubleAreas[t] > limit)) {
            out.set(t);
        }
    }
    return out;
}

// The smoothing rounds with their per point state.
class Smoothing {
public:
    Smoothing(Mesh& mesh, Denoiser::Mode mode, const BitArray& protectedMask, std::span<const int> chainPoints)
        : mesh_(mesh), rules_(rulesFor(mode)), adjacency_(mesh.adjacency()), n_(mesh.pointCount()),
          fixed_(protectedMask.size() == static_cast<std::size_t>(n_) ? protectedMask
                                                                      : BitArray(static_cast<std::size_t>(n_))),
          flat_(static_cast<std::size_t>(n_)), kind_(static_cast<std::size_t>(n_), 0.0f),
          blend_(static_cast<std::size_t>(n_), 0.0f)
    {
        for (const int v : chainPoints) {
            if (v < 0) {
                continue;
            }
            if (static_cast<std::size_t>(v) >= fixed_.size()) {
                fixed_.resize(static_cast<std::size_t>(v) + 1);
            }
            fixed_.set(static_cast<std::size_t>(v));
        }
    }

    void run()
    {
        for (int r = 0; r < rules_.rounds; ++r) {
            // Round 0 reads the edge angles and normals the chain detection left in the mesh.
            // The defects of the recheck stay valid when no point moved since then.
            const std::vector<float>& angles = mesh_.dihedralAngles();
            const std::vector<float> defects =
                keptDefects_ ? std::move(*keptDefects_) : angleDefects(adjacency_, cornerAngles(mesh_));
            keptDefects_.reset();
            const std::vector<float> lowest = neighbourMinimum(defects);
            std::vector<float> weights(static_cast<std::size_t>(n_), 0.0f);
            bool any = false;
            for (int v = 0; v < n_; ++v) {
                any = weigh(v, angles, defects, lowest, weights) || any;
            }
            if (!any) {
                break;
            }
            if (rules_.detection) {
                weights = maxMeanWeights(weights, 0.35f, 0.4f);
                weights = smoothWeights(weights, 0.6f);
            } else {
                weights = smoothWeights(weights, 0.4f);
                weights = smoothWeights(weights, 0.4f);
            }
            for (float& w : weights) {
                if (w > rules_.weightCap) {
                    w = rules_.weightCap;
                }
            }
            for (const int v : adjacency_.boundaryPoints()) {
                weights[static_cast<std::size_t>(v)] = 0.0f;
            }
            if (rules_.zeroNonManifold) {
                for (const int e : adjacency_.nonManifoldEdges()) {
                    for (const int v : adjacency_.edgeVertices(e)) {
                        weights[static_cast<std::size_t>(v)] = 0.0f;
                    }
                }
            }

            const TriangleGeometry& geometry = mesh_.triangleGeometry();
            Smoother::Options options;
            options.iterations = (r == 0 && !rules_.detection) ? 2 : 3;
            options.foldGuard = true;
            options.guardPoints = &flat_;
            std::vector<Vec3f> positions = Smoother(mesh_, options).applyAll(weights, geometry.normals());
            const std::vector<Vec3f> oldNormals = geometry.normals();
            const BitArray degenerate = degenerateTriangles(geometry);
            mesh_.swapPoints(positions);
            const std::vector<Vec3f>& before = positions;
            if (recheck(oldNormals, degenerate)) {
                revert(before, oldNormals, degenerate);
            }
        }
    }

private:
    [[nodiscard]] std::vector<float> neighbourMinimum(const std::vector<float>& defects) const
    {
        std::vector<float> out(defects.size());
        for (int v = 0; v < n_; ++v) {
            float m = defects[static_cast<std::size_t>(v)];
            for (const int u : adjacency_.pointNeighbours(v)) {
                const float d = defects[static_cast<std::size_t>(u)];
                m = d < m ? d : m;
            }
            out[static_cast<std::size_t>(v)] = m;
        }
        return out;
    }

    void fix(int v, std::vector<float>& weights)
    {
        weights[static_cast<std::size_t>(v)] = 0.0f;
        fixed_.set(static_cast<std::size_t>(v));
    }

    // Weight of one point; true when it got a weight that counts as movement.
    bool weigh(int v, const std::vector<float>& angles, const std::vector<float>& defects,
               const std::vector<float>& lowest, std::vector<float>& weights)
    {
        const auto idx = static_cast<std::size_t>(v);
        float& weight = weights[idx];
        if (fixed_.testBounded(static_cast<std::size_t>(v))) {
            weight = 0.0f;
            return false;
        }
        const float d = defects[idx];
        const float a = std::fabs(d);
        float pos = 0.0f;
        float neg = 0.0f;
        float angleSum = 0.0f;
        int above178 = 0, above60 = 0, above40 = 0, above20 = 0;
        std::vector<int> steep;
        const auto edges = adjacency_.pointEdges(v);
        const auto neighbours = adjacency_.pointNeighbours(v);
        for (std::size_t k = 0; k < edges.size(); ++k) {
            const auto u = static_cast<std::size_t>(neighbours[k]);
            if (defects[u] > 0.0f) {
                pos = pos + defects[u];
            } else {
                neg = neg - lowest[u];
            }
            const float t = angles[static_cast<std::size_t>(edges[k])];
            angleSum = angleSum + t;
            if (wide(t) > degD(178)) {
                ++above178;
            }
            if (wide(t) > degD(60)) {
                ++above60;
                steep.push_back(edges[k]);
            } else if (wide(t) > degD(40)) {
                ++above40;
                steep.push_back(edges[k]);
            } else if (wide(t) > degD(20)) {
                ++above20;
            }
        }
        // Flatness class in float; the weight rules below use the 20 degree limit in double.
        if (a < degF(20)) {
            flat_.set(idx);
            if (a < degF(10)) {
                blend_[idx] = 0.5f;
                kind_[idx] = 0.0f;
            } else {
                kind_[idx] = 0.5f;
            }
        } else {
            kind_[idx] = 1.0f;
        }
        if (above178 > 0) {
            weight = 1.0f;
            return true;
        }
        if (wide(a) < degD(2) || wide(neg) < degD(5) || wide(pos) < degD(5)) {
            if (wide(angleSum) < degD(2)) {
                if (rules_.detection) {
                    weight = 0.5f;
                    return true;
                }
                fix(v, weights);
                return false;
            }
            const bool uneven = above20 >= 2 || (d > 0.0f && wide(pos) < degD(1) && wide(neg) > degD(5)) ||
                                (d < 0.0f && wide(pos) > degD(5) && wide(neg) < degD(1)) ||
                                above20 + above40 + above60 > 4;
            if (uneven) {
                weight = 0.3f;
                return true;
            }
            if (rules_.detection && !(above60 == 2 && wide(a) < degD(4))) {
                weight = 0.0f;
                return false;
            }
            fix(v, weights);
            return false;
        }
        float g;
        if (d > 0.0f) {
            if (neg > (pos + d) * 0.8f) {
                g = ramp(d);
            } else {
                const float q = neg / a;
                const float m = 1.0f <= q ? 1.0f : q;
                g = ramp(m * d);
                if (wide(g) < 0.005) {
                    weight = 0.0f;
                    return false;
                }
            }
        } else {
            if (wide(pos) < degD(0.3)) {
                weight = 0.0f;
                return false;
            }
            if (pos > (neg - d) * 0.8f) {
                g = 1.0f;
            } else {
                const float q = pos / a;
                const float m = 1.0f <= q ? 1.0f : q;
                g = ramp(m * d);
                if (wide(g) < 0.005) {
                    weight = 0.0f;
                    return false;
                }
            }
        }
        // Cube corner: three edges near 90 degrees meet at a defect near 90 degrees.
        auto nearRight = [](float x) { return wide(std::fabs(x - degF(90))) < degD(10); };
        if (nearRight(d) && steep.size() == 3 && above40 == 0) {
            bool corner = true;
            for (const int e : steep) {
                corner = corner && nearRight(angles[static_cast<std::size_t>(e)]);
            }
            if (corner) {
                fix(v, weights);
                return false;
            }
        }
        if (d > degF(90) || !(wide(a) < degD(20))) {
            weight = g;
            return true;
        }
        if (rules_.detection) {
            weight = (1.0f * 0.6f) * g;
            return true;
        }
        fix(v, weights);
        return false;
    }

    // Jacobi blend of every free point's weight with the mean of its neighbours (share f).
    [[nodiscard]] std::vector<float> smoothWeights(const std::vector<float>& in, float f) const
    {
        std::vector<float> out(in.size());
        for (int v = 0; v < n_; ++v) {
            const auto idx = static_cast<std::size_t>(v);
            const auto neighbours = adjacency_.pointNeighbours(v);
            if (fixed_.testBounded(static_cast<std::size_t>(v)) || neighbours.empty()) {
                out[idx] = in[idx];
                continue;
            }
            float sum = 0.0f;
            for (const int u : neighbours) {
                sum += in[static_cast<std::size_t>(u)];
            }
            out[idx] = in[idx] * (1.0f - f) + (sum / static_cast<float>(neighbours.size())) * f;
        }
        return out;
    }

    // Jacobi blend of every free point's weight with the largest weight of it
    // and its neighbours (share a) and with their mean (share b).
    [[nodiscard]] std::vector<float> maxMeanWeights(const std::vector<float>& in, float a, float b) const
    {
        std::vector<float> out(in.size());
        for (int v = 0; v < n_; ++v) {
            const auto idx = static_cast<std::size_t>(v);
            const auto neighbours = adjacency_.pointNeighbours(v);
            if (fixed_.testBounded(static_cast<std::size_t>(v)) || neighbours.empty()) {
                out[idx] = in[idx];
                continue;
            }
            float sum = 0.0f;
            float largest = in[idx];
            for (const int u : neighbours) {
                const float w = in[static_cast<std::size_t>(u)];
                sum += w;
                largest = std::max(largest, w);
            }
            out[idx] = in[idx] * ((1.0f - b) - a) + largest * a + (sum / static_cast<float>(neighbours.size())) * b;
        }
        return out;
    }

    [[nodiscard]] bool flipped(const std::vector<Vec3f>& normals, const std::vector<Vec3f>& oldNormals,
                               const BitArray& degenerate, int t) const
    {
        if (degenerate.testBounded(static_cast<std::size_t>(t))) {
            return false;
        }
        const auto i = static_cast<std::size_t>(t);
        return flipLimit > dot(normals[i], oldNormals[i]);
    }

    // Flipped triangles: one free, not quite flat point each (the one with the
    // largest new angle defect) is smoothed once more. True when that happened.
    bool recheck(const std::vector<Vec3f>& oldNormals, const BitArray& degenerate)
    {
        const std::vector<Vec3f>& normals = mesh_.triangleGeometry().normals();
        const int triangles = mesh_.faceCount();
        std::vector<float> weights(static_cast<std::size_t>(n_), 0.0f);
        int count = 0;
        std::optional<std::vector<float>> defects;
        constexpr float none = -1000.0f;
        for (int t = 0; t < triangles; ++t) {
            if (!flipped(normals, oldNormals, degenerate, t)) {
                continue;
            }
            if (!defects) {
                // The points just moved: no corner angles are known, the defects come from
                // the angles at the corners in double.
                defects = directAngleDefects(mesh_);
            }
            float best = none;
            int chosen = -1;
            for (int i = 0; i < 3; ++i) {
                const int u = adjacency_.vertex(3 * t + i);
                const float du = (*defects)[static_cast<std::size_t>(u)];
                if (!fixed_.testBounded(static_cast<std::size_t>(u)) && kind_[static_cast<std::size_t>(u)] != 0.0f && du > best) {
                    best = du;
                    chosen = u;
                }
            }
            if (chosen != -1 && weights[static_cast<std::size_t>(chosen)] == 0.0f) {
                weights[static_cast<std::size_t>(chosen)] = 0.8f;
                ++count;
            }
        }
        if (count == 0) {
            keptDefects_ = std::move(defects);
            return false;
        }
        Smoother::Options options;
        options.iterations = 3;
        std::vector<Vec3f> positions = Smoother(mesh_, options).applyAll(weights, {});
        mesh_.swapPoints(positions);
        return true;
    }

    // Up to three passes over the still flipped triangles: their points with a
    // blend factor go halfway back, completely back in the second pass.
    void revert(const std::vector<Vec3f>& before, const std::vector<Vec3f>& oldNormals, const BitArray& degenerate)
    {
        for (int pass = 1; pass <= 3; ++pass) {
            const std::vector<Vec3f>& normals = mesh_.triangleGeometry().normals();
            const int triangles = mesh_.faceCount();
            BitArray listed(static_cast<std::size_t>(n_));
            std::vector<int> list;
            bool any = false;
            for (int t = 0; t < triangles; ++t) {
                if (!flipped(normals, oldNormals, degenerate, t)) {
                    continue;
                }
                any = true;
                for (int i = 0; i < 3; ++i) {
                    const int u = adjacency_.vertex(3 * t + i);
                    const auto idx = static_cast<std::size_t>(u);
                    if (!fixed_.testBounded(static_cast<std::size_t>(u)) && blend_[idx] != 0.0f && !listed.test(idx)) {
                        list.push_back(u);
                        listed.set(idx);
                    }
                }
            }
            if (!any) {
                return;
            }
            const std::vector<Vec3f>& current = mesh_.points();
            std::vector<Vec3f> moved;
            moved.reserve(list.size());
            for (const int u : list) {
                const auto idx = static_cast<std::size_t>(u);
                const Vec3f old = before[idx];
                if (pass == 2) {
                    moved.push_back(old);
                } else {
                    const Vec3f p = current[idx];
                    const float z = blend_[idx];
                    const float keep = 1.0f - z;
                    moved.push_back({p.x * keep + old.x * z, p.y * keep + old.y * z, p.z * keep + old.z * z});
                }
            }
            mesh_.movePoints(list, moved);
        }
    }

    Mesh& mesh_;
    ModeRules rules_;
    const Adjacency& adjacency_;
    int n_;
    BitArray fixed_; // weight 0 for good, grows over the rounds
    BitArray flat_;  // angle defect below 20 degrees in some round; the fold guard works on these
    std::vector<float> kind_;  // 1, 0.5 or 0 by angle defect; 0 keeps a point out of the recheck
    std::vector<float> blend_; // 0.5 once the defect was below 10 degrees: share of the old position when reverting
    std::optional<std::vector<float>> keptDefects_; // defects of a recheck that moved nothing, for the next round
};

} // namespace

FeatureChains::Result Denoiser::run(Mesh& mesh, const Options& options)
{
    // Chain detection on normals with replaced degenerate triangles; switching
    // the mode back keeps the replaced normals (and the edge angles from them)
    // until the points move.
    (void)mesh.triangleGeometry();
    const int previousMode = mesh.replacementNormalMode();
    mesh.setReplacementNormalMode(2, true);
    FeatureChains::Result chains = FeatureChains::classify(mesh, options.chains);
    mesh.setReplacementNormalMode(previousMode, false);

    const BitArray none;
    Smoothing(mesh, options.mode, options.fixedPoints != nullptr ? *options.fixedPoints : none, chains.points.entries)
        .run();
    return chains;
}

void Denoiser::run(Mesh& mesh, const PrepState& state)
{
    if (mesh.topology().mode() != FaceMode::Triangles) {
        return;
    }
    const BitArray protectedMask = protectedPoints(mesh);
    Options options;
    options.chains.creaseAngle = state.creaseAngle;
    // Without a stored area the area of the mesh (the reference lengths would be NaN otherwise).
    options.chains.lengths = FeatureChains::referenceLengths(state.knownArea(mesh), state.quadCount);
    options.fixedPoints = &protectedMask;
    (void)run(mesh, options);
}

} // namespace areko::remesher::prep
