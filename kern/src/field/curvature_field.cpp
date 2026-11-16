#include "curvature_field.h"

#include "curvature_radii.h"
#include "curve_influence.h"
#include "edge_curvature.h"
#include "symmetric_eigen3.h"

#include "mesh/measures.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <optional>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;
using curves::GuideCurve;
namespace flag = curvature_flags;

namespace {

constexpr float fromBits(std::uint32_t bits) { return std::bit_cast<float>(bits); }

// Anisotropy a triangle needs for a reliable direction.
constexpr float anisotropyLimit = 0.7f;
static_assert(std::bit_cast<std::uint32_t>(anisotropyLimit) == 0x3f333333u);
// Sines of 30, 25, 42.5 and 60 degrees as sinf gives them for the float angles.
constexpr float sinTangent = fromBits(0x3f000000u);
constexpr float sinLow = fromBits(0x3ed8616cu);
constexpr float sinMid = fromBits(0x3f2cf37bu);
constexpr float sinHigh = fromBits(0x3f5db3d8u);
// Angles of the level choice: 4 and 8 degrees, pi / 2 (float and its double).
constexpr float stableSpread = fromBits(0x3d8efa35u);
constexpr float fallbackSpread = fromBits(0x3e0efa35u);
constexpr float halfPiF = fromBits(0x3fc90fdbu);
constexpr double halfPiOfFloat = 1.5707963705062866;
constexpr double angle5 = 0.08726646259971647;
constexpr double angle20 = 0.3490658503988659;
constexpr float tinyWeight = fromBits(0x358637bdu);
constexpr float fixedWeight = 50000.0f;
static_assert(std::bit_cast<std::uint32_t>(0.39999998f) == 0x3eccccccu);
static_assert(std::bit_cast<std::uint32_t>(59.9f) == 0x426f999au && std::bit_cast<std::uint32_t>(-59.9f) == 0xc26f999au);
static_assert(std::bit_cast<std::uint32_t>(1.5707964f) == 0x3fc90fdbu);

// Results of one measuring level per triangle.
struct Level {
    float radius = 0.0f;
    float minCurvature = 0.0f;
    std::vector<Vec3f> direction;
    std::vector<float> lambda1;
    std::vector<float> lambda2; // with its sign
    BitArray valid;

    // (|l2| - |l1|) / |l2|, 0 for |l2| = 0.
    [[nodiscard]] float anisotropy(std::size_t t) const
    {
        const float k2 = std::fabs(lambda2[t]);
        return k2 == 0.0f ? 0.0f : (k2 - std::fabs(lambda1[t])) / k2;
    }
};

using Tensor = std::array<double, 6>;

// The working state of one curvature computation.
class Computation {
public:
    Computation(parts::Part& part, const CurvatureField::Options& options, const run::Progress& progress)
        : part_(part), mesh_(part.mesh), adjacency_(part.mesh.adjacency()),
          normals_(part.mesh.triangleGeometry().normals()), options_(options), progress_(progress),
          faceCount_(static_cast<std::size_t>(part.mesh.faceCount()))
    {
    }

    void run();

private:
    // Mirror pairs the curvature may use: those of a symmetric part only.
    [[nodiscard]] prep::SymmetryCache* mirrorPairs();
    void markMirrorHalf();
    [[nodiscard]] bool countsInMean(int v) const;
    void fixCurveDirections();
    void fixBorderDirections();
    void fixPlaneDirections();
    void setFixed(int f, Vec3f direction, std::uint8_t mark);
    [[nodiscard]] std::vector<Tensor> pointTensors(float radius) const;
    void triangleDirections(Level& level, const std::vector<Tensor>& tensors) const;
    void selectLevels();
    [[nodiscard]] bool stable(int l, std::size_t t) const;
    [[nodiscard]] bool selectStable(std::size_t t);
    [[nodiscard]] bool selectFixed(std::size_t t);
    [[nodiscard]] bool selectFallback(std::size_t t);
    void averageLevels(std::size_t t);
    void mirrorResults();

    parts::Part& part_;
    const Mesh& mesh_;
    const Adjacency& adjacency_;
    const std::vector<Vec3f>& normals_;
    CurvatureField::Options options_;
    const run::Progress& progress_;
    std::size_t faceCount_;
    CurvatureRadii radii_;
    std::optional<EdgeCurvature> edges_;
    std::vector<int> faceMirror_;
    std::vector<Level> levels_;
    Curvature result_;
};

void Computation::run()
{
    const float area = mesh_.triangleGeometry().totalArea();
    radii_ = CurvatureRadii::compute(part_.edgeLength, area, part_.bordersFrozen, options_.severalParts);
    CurveInfluence::assignWidths(part_.curves, radii_.maxRadius());
    if (radii_.maxRadius() == 0.0f) {
        return;
    }
    markMirrorHalf();
    progress_.check();
    edges_.emplace(mesh_);
    result_.sharpPoints = edges_->sharpPoints(mesh_, part_.curves);
    progress_.check();

    result_.direction.assign(faceCount_, Vec3f{});
    result_.lambda1.assign(faceCount_, 0.0f);
    result_.lambda2.assign(faceCount_, 0.0f);
    result_.flags.assign(faceCount_, 0);
    result_.weight.assign(faceCount_, 1.0f);
    CurveInfluence::prepareBands(part_.curves, mesh_);
    result_.curveDistance = CurveInfluence::distance(mesh_, part_.curves);
    fixCurveDirections();
    progress_.check();
    if (options_.followBorders) {
        fixBorderDirections();
        progress_.check();
    }
    if (options_.symmetry) {
        fixPlaneDirections();
        progress_.check();
    }

    levels_.resize(static_cast<std::size_t>(radii_.levelCount()));
    for (int l = 0; l < radii_.levelCount(); ++l) {
        Level& level = levels_[static_cast<std::size_t>(l)];
        level.radius = radii_.levelRadius(l);
        level.minCurvature = CurvatureRadii::minCurvature(level.radius);
        triangleDirections(level, pointTensors(level.radius));
        progress_.check();
    }
    selectLevels();
    mirrorResults();
    part_.curvature = std::move(result_);
}

prep::SymmetryCache* Computation::mirrorPairs()
{
    return part_.mirror.symmetric(mesh_, true);
}

void Computation::markMirrorHalf()
{
    result_.mirrorFaces = BitArray(faceCount_);
    if (!options_.halfComputed) {
        return;
    }
    prep::SymmetryCache* cache = part_.mirror.evaluate(mesh_, true);
    if (cache == nullptr) {
        return;
    }
    if (!part_.mirror.classes.valid || !cache->allFacesMirrored(mesh_)) {
        return;
    }
    faceMirror_ = cache->faceMirror(mesh_);
    const Topology& topology = mesh_.topology();
    for (std::size_t f = 0; f < faceCount_; ++f) {
        // Only triangles with a mirror triangle are left to the mirror.
        if (cache->sideOfFace(topology, static_cast<int>(f)) == 0 && f < faceMirror_.size() && faceMirror_[f] != -1) {
            result_.mirrorFaces.set(f);
            ++result_.mirrorFaceCount;
        }
    }
    const std::vector<int>& pairs = cache->pointMirror();
    const BitArray& sides = cache->pointSide();
    result_.mirrorPoints = BitArray(static_cast<std::size_t>(mesh_.pointCount()));
    for (std::size_t v = 0; v < pairs.size() && v < result_.mirrorPoints.size(); ++v) {
        const int m = pairs[v];
        if (!sides.testBounded(v) && m != -1 && m != static_cast<int>(v)) {
            result_.mirrorPoints.set(v);
        }
    }
}

bool Computation::countsInMean(int v) const
{
    // Points whose tensor is not computed do not count.
    const std::size_t i = static_cast<std::size_t>(v);
    return !result_.sharpPoints.testBounded(i) && !result_.mirrorPoints.testBounded(i);
}

void Computation::setFixed(int f, Vec3f direction, std::uint8_t mark)
{
    result_.direction[static_cast<std::size_t>(f)] = direction;
    result_.flags[static_cast<std::size_t>(f)] = mark;
}

void Computation::fixCurveDirections()
{
    const std::vector<GuideCurve>& curves = part_.curves;
    const CurveDistance& distance = result_.curveDistance;
    std::vector<std::vector<Vec3f>> tangents(curves.size());
    for (std::size_t i = 0; i < curves.size(); ++i) {
        if (curves[i].active && curves::isBand(curves[i].kind)) {
            tangents[i] = CurveInfluence::tangents(curves[i].cutPoints);
        }
    }
    auto seedTangent = [&](int seed) {
        const std::size_t s = static_cast<std::size_t>(seed);
        return tangents[static_cast<std::size_t>(distance.seedCurve[s])][static_cast<std::size_t>(distance.seedPoint[s])];
    };
    const Topology& topology = mesh_.topology();
    for (std::size_t t = 0; t < faceCount_; ++t) {
        if (distance.faceDistance[t] > 0.5f) {
            continue;
        }
        int corner[3];
        bool bands = true;
        for (int k = 0; k < 3 && bands; ++k) {
            corner[k] = topology.vertex(3 * static_cast<int>(t) + k);
            const int seed = distance.pointSeed[static_cast<std::size_t>(corner[k])];
            bands = seed >= 0 && curves::isBand(curves[static_cast<std::size_t>(distance.seedCurve[static_cast<std::size_t>(seed)])].kind);
        }
        if (!bands) {
            continue;
        }
        int nearest = corner[0];
        for (int k = 1; k < 3; ++k) {
            if (distance.pointDistance[static_cast<std::size_t>(corner[k])] < distance.pointDistance[static_cast<std::size_t>(nearest)]) {
                nearest = corner[k];
            }
        }
        const Vec3f tangent = seedTangent(distance.pointSeed[static_cast<std::size_t>(nearest)]);
        setFixed(static_cast<int>(t), projectOntoPlane(tangent, normals_[t]), flag::fixedCurve);
    }
    const std::vector<Vec3f>& points = mesh_.points();
    for (const GuideCurve& curve : curves) {
        if (!curve.active || !curves::cutsIntoMesh(curve.kind)) {
            continue;
        }
        const std::vector<int>& chain = curve.cutVertices;
        for (std::size_t j = 0; j + 1 < chain.size(); ++j) {
            const int a = chain[j];
            const int b = chain[j + 1];
            const int e = adjacency_.edgeBetween(a, b);
            if (e == -1) {
                continue;
            }
            const auto corners = adjacency_.edgeCorners(e);
            if (corners.size() > 2) {
                continue;
            }
            Vec3f d = points[static_cast<std::size_t>(b)] - points[static_cast<std::size_t>(a)];
            normalizeInPlace(d);
            for (const int c : corners) {
                setFixed(c / 3, d, flag::fixedCurve);
            }
        }
    }
}

void Computation::fixBorderDirections()
{
    for (const int e : adjacency_.boundaryEdges()) {
        for (const int c : adjacency_.edgeCorners(e)) {
            setFixed(c / 3, edges_->directions()[static_cast<std::size_t>(e)], flag::fixedBorder);
        }
    }
}

void Computation::fixPlaneDirections()
{
    const prep::SymmetryCache* cache = mirrorPairs();
    if (cache == nullptr || cache->pointMirror().empty()) {
        return;
    }
    const std::vector<int>& pairs = cache->pointMirror();
    auto onPlane = [&pairs](int v) {
        return static_cast<std::size_t>(v) < pairs.size() && pairs[static_cast<std::size_t>(v)] == v;
    };
    for (int e = 0; e < adjacency_.edgeCount(); ++e) {
        const auto ends = adjacency_.edgeVertices(e);
        if (!onPlane(ends[0]) || !onPlane(ends[1])) {
            continue;
        }
        for (const int c : adjacency_.edgeCorners(e)) {
            setFixed(c / 3, edges_->directions()[static_cast<std::size_t>(e)], flag::fixedPlane);
        }
    }
}

std::vector<Tensor> Computation::pointTensors(float radius) const
{
    const std::size_t pointCount = static_cast<std::size_t>(mesh_.pointCount());
    std::vector<Tensor> out(pointCount, Tensor{});
    if (!(0.0f < radius)) {
        return out;
    }
    const std::vector<Vec3f>& points = mesh_.points();
    const std::vector<float>& lengths = edges_->lengths();
    const std::vector<float>& angles = edges_->angles();
    const BitArray& sharp = result_.sharpPoints;
    const float r2 = radius * radius;
    const double scale = static_cast<double>(1.0f / ((radius * radius) * piF));
    std::vector<std::uint32_t> pointMark(pointCount, 0);
    std::vector<std::uint32_t> edgeMark(static_cast<std::size_t>(adjacency_.edgeCount()), 0);
    std::uint32_t mark = 0;
    std::vector<int> stack;
    for (std::size_t v = 0; v < pointCount; ++v) {
        if (sharp.testBounded(v) || result_.mirrorPoints.testBounded(v)) {
            continue;
        }
        const Vec3f p = points[v];
        ++mark;
        pointMark[v] = mark;
        stack.assign(1, static_cast<int>(v));
        Tensor sum{};
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            const auto pointEdges = adjacency_.pointEdges(u);
            const auto neighbours = adjacency_.pointNeighbours(u);
            for (std::size_t k = 0; k < pointEdges.size(); ++k) {
                const int e = pointEdges[k];
                const std::size_t ei = static_cast<std::size_t>(e);
                if (edgeMark[ei] == mark) {
                    continue;
                }
                edgeMark[ei] = mark;
                if (lengths[ei] == 0.0f) {
                    continue;
                }
                const int w = neighbours[k];
                const Vec3f dw = points[static_cast<std::size_t>(w)] - p;
                const float d2 = dot(dw, dw);
                const std::span<const double, 6> tensor = edges_->tensor(e);
                if (r2 <= d2) {
                    // The edge leaves the sphere: its share by the radial distances.
                    if (angles[ei] != 0.0f) {
                        const Vec3f du = points[static_cast<std::size_t>(u)] - p;
                        const float outer = std::sqrt(d2);
                        const float inner = std::sqrt(dot(du, du));
                        const double share = static_cast<double>((radius - inner) / (outer - inner));
                        for (std::size_t i = 0; i < 6; ++i) {
                            sum[i] = sum[i] + tensor[i] * share;
                        }
                    }
                    continue;
                }
                const std::size_t wi = static_cast<std::size_t>(w);
                if (pointMark[wi] != mark) {
                    pointMark[wi] = mark;
                    if (!sharp.testBounded(wi)) {
                        stack.push_back(w);
                    }
                }
                if (angles[ei] != 0.0f) {
                    for (std::size_t i = 0; i < 6; ++i) {
                        sum[i] = sum[i] + tensor[i];
                    }
                }
            }
        }
        for (std::size_t i = 0; i < 6; ++i) {
            out[v][i] = sum[i] * scale;
        }
    }
    return out;
}

void Computation::triangleDirections(Level& level, const std::vector<Tensor>& tensors) const
{
    level.direction.assign(faceCount_, Vec3f{});
    level.lambda1.assign(faceCount_, 0.0f);
    level.lambda2.assign(faceCount_, 0.0f);
    level.valid = BitArray(faceCount_);
    const Topology& topology = mesh_.topology();
    const std::vector<float>& curveDistance = result_.curveDistance.faceDistance;
    for (std::size_t f = 0; f < faceCount_; ++f) {
        if (result_.mirrorFaces.testBounded(f)) {
            continue;
        }
        Tensor mean{};
        int n = 0;
        for (int c = 3 * static_cast<int>(f); c < 3 * static_cast<int>(f) + 3; ++c) {
            const int v = topology.vertex(c);
            if (!countsInMean(v)) {
                continue;
            }
            const Tensor& t = tensors[static_cast<std::size_t>(v)];
            for (std::size_t i = 0; i < 6; ++i) {
                mean[i] = n == 0 ? t[i] : mean[i] + t[i];
            }
            ++n;
        }
        if (n >= 1) {
            const double s = 1.0 / static_cast<double>(n);
            for (double& m : mean) {
                m = m * s;
            }
        }
        const std::array<SymmetricEigen3::Pair, 3> pairs = SymmetricEigen3::solve(mean);
        const Vec3d& v2 = pairs[2].vector;
        const Vec3f d{static_cast<float>(v2.x), static_cast<float>(v2.y), static_cast<float>(v2.z)};
        level.direction[f] = d;
        level.lambda1[f] = static_cast<float>(pairs[1].value);
        level.lambda2[f] = static_cast<float>(pairs[2].value);
        const float k2 = std::fabs(level.lambda2[f]);
        float limit = anisotropyLimit;
        const float dn = curveDistance[f];
        if (0.5f <= dn && dn < 1.5f) {
            // Near a curve less anisotropy is needed.
            limit = (limit - 1.0f) * (dn - 0.5f) + 1.0f;
        }
        const float a = level.anisotropy(f);
        const float tangential = std::fabs(dot(normals_[f], d));
        if (!(limit > a) && !(level.minCurvature > k2) && !(tangential > sinTangent)) {
            level.valid.set(f);
        }
    }
}

bool Computation::stable(int l, std::size_t t) const
{
    const int lo = std::max(0, l - 1);
    const int hi = std::min(l + 1, radii_.levelCount() - 1);
    for (int j = lo; j <= hi; ++j) {
        if (!levels_[static_cast<std::size_t>(j)].valid.test(t)) {
            return false;
        }
    }
    return true;
}

void Computation::selectLevels()
{
    for (std::size_t t = 0; t < faceCount_; ++t) {
        if (result_.mirrorFaces.testBounded(t)) {
            continue;
        }
        const bool fixed = (result_.flags[t] & flag::fixedMask) != 0;
        if (fixed ? selectFixed(t) : selectStable(t)) {
            continue;
        }
        if (!selectFallback(t)) {
            averageLevels(t);
        }
    }
}

// A free triangle takes the stable level whose direction varies least
// against its neighbouring levels.
bool Computation::selectStable(std::size_t t)
{
    const int n = radii_.levelCount();
    int best = -1;
    float bestSpread = FLT_MAX;
    for (int l = 0; l < n; ++l) {
        if (!levels_[static_cast<std::size_t>(l)].valid.test(t) || !stable(l, t)) {
            continue;
        }
        const int lo = std::max(0, l - 1);
        const int hi = std::min(l + 1, n - 1);
        const Vec3f d = levels_[static_cast<std::size_t>(l)].direction[t];
        float sum = 0.0f;
        for (int j = lo; j <= hi; ++j) {
            if (j != l) {
                sum = sum + lineAngle(d, levels_[static_cast<std::size_t>(j)].direction[t]);
            }
        }
        const float spread = sum / static_cast<float>(hi - lo + 1);
        if (spread < bestSpread) {
            best = l;
            bestSpread = spread;
        }
    }
    if (best == -1) {
        return false;
    }
    const Level& level = levels_[static_cast<std::size_t>(best)];
    result_.direction[t] = level.direction[t];
    result_.lambda1[t] = level.lambda1[t];
    result_.lambda2[t] = level.lambda2[t];
    result_.flags[t] = static_cast<std::uint8_t>(result_.flags[t] | flag::direction | flag::reliable);
    const float ws = bestSpread > stableSpread ? 500.0f : (-500.0f * bestSpread) / stableSpread + 1000.0f;
    const float a = level.anisotropy(t);
    float wa;
    if (a < anisotropyLimit) {
        wa = 1500.0f;
    } else if (a > 1.0f) {
        wa = 3000.0f;
    } else {
        wa = (((a - anisotropyLimit) * 500.0f) / (1.0f - anisotropyLimit) + 500.0f) * 3.0f;
    }
    result_.weight[t] = (wa + ws) * 0.25f;
    return true;
}

// A triangle with a fixed direction keeps it and takes the curvature values
// of the stable level that agrees best with it, blended by their agreement.
bool Computation::selectFixed(std::size_t t)
{
    result_.weight[t] = fixedWeight;
    const int n = radii_.levelCount();
    const Vec3f fixedDirection = result_.direction[t];
    int best = -1;
    float bestCost = FLT_MAX;
    for (int l = 0; l < n; ++l) {
        if (!stable(l, t)) {
            continue;
        }
        const Level& level = levels_[static_cast<std::size_t>(l)];
        const float phi = lineAngle(fixedDirection, level.direction[t]);
        float c = 0.0f;
        if (wide(phi) > angle5) {
            c = phi / halfPiF;
            if (wide(phi) > angle20) {
                c = static_cast<float>(wide(c) + ((wide(phi) - angle20) * 3.0) / halfPiOfFloat);
            }
        }
        const float a = level.anisotropy(t);
        float cost = (a + a) + c;
        if (l == 0) {
            cost = cost * 0.5f;
        } else if (l > n / 2) {
            cost = l == n - 1 ? cost / 1.5f : cost * 1.5f;
        }
        if (cost < bestCost) {
            best = l;
            bestCost = cost;
        }
    }
    if (best == -1) {
        return false;
    }
    const Level& level = levels_[static_cast<std::size_t>(best)];
    float l1 = level.lambda1[t];
    float l2 = level.lambda2[t];
    float c = std::fabs(dot(fixedDirection, level.direction[t]));
    const float mid = (l1 + l2) * 0.5f;
    if (c < 0.5f) {
        std::swap(l1, l2);
        c = 1.0f - c;
    }
    const float u = c - 1.0f;
    result_.lambda1[t] = ((mid - l1) * u) * -2.0f + l1;
    result_.lambda2[t] = ((mid - l2) * u) * -2.0f + l2;
    result_.flags[t] = static_cast<std::uint8_t>(result_.flags[t] | flag::reliable);
    return true;
}

// Without a stable level the cheapest level by anisotropy, curvature,
// spread and tangential share is taken, small radii favoured.
bool Computation::selectFallback(std::size_t t)
{
    const int n = radii_.levelCount();
    std::vector<float> aniso(static_cast<std::size_t>(n));
    for (int l = 0; l < n; ++l) {
        aniso[static_cast<std::size_t>(l)] = levels_[static_cast<std::size_t>(l)].anisotropy(t);
    }
    const Vec3f normal = normals_[t];
    int best = -1;
    float bestCost = FLT_MAX;
    float chosenCost = 0.0f;
    for (int l = 0; l < n; ++l) {
        const Level& level = levels_[static_cast<std::size_t>(l)];
        const float lambda = level.lambda2[t];
        const float k = level.minCurvature;
        const float a = aniso[static_cast<std::size_t>(l)];
        float fa;
        if (a < 0.3f) {
            fa = 1.0f;
        } else if (a > 0.7f) {
            fa = 0.0f;
        } else {
            fa = -(a - 0.3f) / 0.39999998f + 1.0f;
        }
        if (a < 0.05f) {
            fa = fa + 2.0f;
        }
        float fk;
        const bool belowFifth = lambda < k / 5.0f;
        if (belowFifth && l != n - 1) {
            continue;
        }
        if (lambda < k * 0.5f) {
            fk = belowFifth ? 4.0f : ((lambda - k / 5.0f) * -2.5f) / (k * 0.5f - k / 5.0f) + 4.0f;
        } else if (lambda > k) {
            fk = 0.0f;
        } else {
            fk = ((lambda - k * 0.5f) * -1.5f) / (k - k * 0.5f) + 1.5f;
        }
        const float s = fk + fa;
        float fd = 1.0f;
        if (s < 1.5f) {
            // Spread of the direction against all other levels. A level that is
            // not anisotropic itself, or has no other level to compare with, has
            // no spread and is no candidate; a zero spread instead would make
            // weak levels cheap (measured: more triangles and self-intersections
            // on the base model).
            if (wide(a) < 0.3) {
                continue;
            }
            float sum = 0.0f;
            float count = 0.0f;
            for (int j = 0; j < n; ++j) {
                if (j == l) {
                    continue;
                }
                sum = sum + lineAngle(level.direction[t], levels_[static_cast<std::size_t>(j)].direction[t]);
                count = count + 1.0f;
            }
            if (!(count > 0.0f)) {
                continue;
            }
            fd = (sum / count) / halfPiF;
        }
        const float x = std::fabs(dot(normal, level.direction[t]));
        float ft;
        if (x > sinHigh) {
            ft = 100.0f;
        } else if (x > sinMid) {
            ft = ((x - sinMid) * 15.0f) / (sinHigh - sinMid) + 5.0f;
        } else if (x > sinLow) {
            ft = ((x - sinLow) * 5.0f) / (sinMid - sinLow);
        } else {
            ft = 0.0f;
        }
        float cost = (s + fd) + ft;
        const float g = -static_cast<float>(l) / static_cast<float>(n - 1) + 1.5f;
        cost = cost * g;
        if (l == n - 2) {
            cost = cost * 0.5f;
        } else if (l == n - 1) {
            cost = cost * 0.25f;
        }
        if (cost < bestCost) {
            best = l;
            bestCost = cost;
            chosenCost = cost;
        }
    }
    if (best == -1) {
        return false;
    }
    const Level& level = levels_[static_cast<std::size_t>(best)];
    const bool fixed = (result_.flags[t] & flag::direction) != 0;
    if (!fixed) {
        result_.direction[t] = level.direction[t];
        result_.flags[t] = static_cast<std::uint8_t>(result_.flags[t] | flag::direction);
    }
    const float lambda = level.lambda2[t];
    const float k = level.minCurvature;
    const float a = aniso[static_cast<std::size_t>(best)];
    result_.lambda1[t] = level.lambda1[t];
    result_.lambda2[t] = lambda;
    if (wide(k) / 1.5 > wide(lambda) && a > 0.5f) {
        result_.flags[t] = static_cast<std::uint8_t>(result_.flags[t] | flag::reliable);
    }
    if (fixed) {
        result_.weight[t] = fixedWeight;
        return true;
    }
    const float x = std::fabs(dot(normal, level.direction[t]));
    if (a < 0.1f || lambda < k / 20.0f || x > sinHigh) {
        result_.weight[t] = tinyWeight;
        return true;
    }
    float wa;
    if (a < 0.3f) {
        wa = tinyWeight;
    } else if (a < anisotropyLimit) {
        wa = ((a - 0.3f) * 0.9f) / (anisotropyLimit - 0.3f) + 0.1f;
    } else if (a > 1.0f) {
        wa = 60.0f;
    } else {
        wa = ((a - anisotropyLimit) * 59.0f) / (1.0f - anisotropyLimit) + 1.0f;
    }
    float wk;
    if (lambda < k / 2.5f) {
        wk = 0.1f;
    } else if (lambda > k) {
        wk = 60.0f;
    } else {
        wk = ((lambda - k / 2.5f) * 59.9f) / (k - k / 2.5f) + 0.1f;
    }
    float ws;
    if (a < 0.3f || lambda < k / 5.0f) {
        ws = tinyWeight;
    } else if (chosenCost < 0.0f) {
        ws = 60.0f;
    } else if (chosenCost > fallbackSpread) {
        ws = 0.1f;
    } else {
        ws = (chosenCost * -59.9f) / fallbackSpread + 60.0f;
    }
    float w = static_cast<float>(std::pow(wide((wa * ws) * wk), 0.3330000042915344));
    if (x > sinMid) {
        w = w / (((x - sinMid) * 990.0f) / (sinHigh - sinMid) + 10.0f);
    } else if (x > sinLow) {
        w = w / (((x - sinLow) * 9.0f) / (sinMid - sinLow) + 1.0f);
    }
    result_.weight[t] = w;
    return true;
}

// No level could be chosen: the curvature values are the weighted mean of
// their magnitudes over all levels; a free triangle gets almost no weight.
void Computation::averageLevels(std::size_t t)
{
    const int n = radii_.levelCount();
    float total = 0.0f;
    float sum1 = 0.0f;
    float sum2 = 0.0f;
    for (int j = 0; j < n; ++j) {
        const Level& level = levels_[static_cast<std::size_t>(j)];
        const bool valid = level.valid.test(t);
        float w = valid ? 5.0f : 1.0f;
        if (j == 0) {
            w = valid ? 2.5f : 0.5f;
        }
        if (j == n - 1) {
            w = w * 0.5f;
        }
        total = total + w;
        sum1 = sum1 + std::fabs(level.lambda1[t]) * w;
        sum2 = sum2 + std::fabs(level.lambda2[t]) * w;
    }
    result_.lambda1[t] = sum1 / total;
    result_.lambda2[t] = sum2 / total;
    if ((result_.flags[t] & flag::fixedMask) == 0) {
        result_.weight[t] = tinyWeight;
    }
}

void Computation::mirrorResults()
{
    if (result_.mirrorFaceCount == 0) {
        return;
    }
    for (std::size_t f = 0; f < faceCount_; ++f) {
        if (!result_.mirrorFaces.test(f)) {
            continue;
        }
        const int g = faceMirror_[f];
        if (g == -1) {
            continue;
        }
        const std::size_t gi = static_cast<std::size_t>(g);
        result_.lambda1[f] = result_.lambda1[gi];
        result_.lambda2[f] = result_.lambda2[gi];
        result_.flags[f] = result_.flags[gi];
        result_.weight[f] = result_.weight[gi];
    }
}

} // namespace

void CurvatureField::compute(parts::Part& part, const Options& options, const run::Progress& progress)
{
    part.curvature.reset();
    Computation(part, options, progress).run();
}

} // namespace areko::remesher::field
