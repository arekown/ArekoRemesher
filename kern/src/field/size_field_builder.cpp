#include "size_field_builder.h"

#include "border_sizes.h"
#include "density_factors.h"
#include "size_smoothing.h"

#include "solver/direct_solver.h"
#include "solver/normal_equations.h"
#include "util/select.h"

#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <optional>
#include <utility>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;
using util::maxsd;
using util::maxss;
using util::minss;

namespace {

// Half a turn: a size is the chord of one of the n segments of a full turn, 2 sin(halfTurn / n) over the curvature.
constexpr double halfTurn = piD;
// Chord of five segments per full turn, 2 sin(pi / 5): no size is pushed below it by the floor.
constexpr float fiveSegments = 1.1755705f;
static_assert(std::bit_cast<std::uint32_t>(fiveSegments) == 0x3f967918u);

// A direction counts as reliable above this trust weight; curvatures below the threshold count as flat.
constexpr float trustLimit = 0.01f;
constexpr double flatLimit = 1e-8;
// Share of the smallest size the sizes are never pushed below.
constexpr float floorShare = 0.9f;
// Band of estimated quads over the budget that ends the counting.
constexpr float bandTop = 1.05f;
constexpr double bandBottom = 0.97;
// Data weights of the fit: base, gain towards the smallest size, weight across a frozen border.
constexpr float dataWeight = 10.0f;
constexpr float dataGain = 100.0f;
constexpr double borderWeight = 1000.0;
// Weight of an edge row: at least edgeWeightMin, else the mean data weight of its two sides times edgeWeightFactor.
constexpr double edgeWeightMin = 20.0;
constexpr double edgeWeightFactor = 6.0;

// The working state of one size field computation.
class Computation {
public:
    Computation(parts::Part& part, const SizeFieldBuilder::Options& options, const run::Progress& progress)
        : part_(part), mesh_(part.mesh), adjacency_(part.mesh.adjacency()), curvature_(*part.curvature),
          cuts_(*part.cutGraph), options_(options), progress_(progress),
          faceCount_(static_cast<std::size_t>(part.mesh.faceCount())), marks_(part.frozenBorderMarks()),
          smoothing_(adjacency_, cuts_)
    {
    }

    void run();

private:
    // The quad budget is kept strictly. Frozen borders decide the sizes
    // themselves: such a part does not count, whatever the setting says.
    [[nodiscard]] bool exactCount() const noexcept { return options_.exactCount && marks_ == nullptr; }
    [[nodiscard]] bool inMirrorHalf(std::size_t t) const noexcept { return cuts_.mirrorHalf.testBounded(t); }

    void setLimits();
    void fill(float size);
    // The sizes from curvature, smoothing and counting; false when the field is final without the fit.
    [[nodiscard]] bool adapt();
    void followCurvature(float pointsPerCircle, float floor);
    void mirror();
    void smooth();
    void clamp(float floor);
    void fit();
    void applyDensity();

    // Rows of the fit: per row its entries (unknown, value), target and weight.
    struct Rows {
        std::vector<solver::ListMatrix::List> entries;
        std::vector<double> targets;
        std::vector<double> weights;
    };
    using Unknowns = std::vector<std::array<int, 2>>;
    using Weights = std::vector<std::array<double, 2>>;
    [[nodiscard]] Unknowns numberUnknowns(int& count) const;
    [[nodiscard]] Weights dataWeights() const;
    [[nodiscard]] float size(std::size_t t, int component) const { return (component == 0 ? field_.sizeU : field_.sizeV)[t]; }
    void addEdgeRow(const Unknowns& unknowns, const Weights& weights, std::size_t t0, int p, std::size_t t1, int q,
                    Rows& rows) const;

    parts::Part& part_;
    const Mesh& mesh_;
    const Adjacency& adjacency_;
    const Curvature& curvature_;
    const CutGraph& cuts_;
    SizeFieldBuilder::Options options_;
    const run::Progress& progress_;
    std::size_t faceCount_;
    // Freeze marks of a part with frozen borders, else null.
    const std::vector<FreezeMark>* marks_;
    SizeSmoothing smoothing_;
    std::optional<BorderSizes> border_;
    SizeField field_;
};

void Computation::run()
{
    field_.edgeLength = part_.edgeLength;
    setLimits();
    progress_.check();
    if (marks_ != nullptr) {
        border_.emplace(part_, *marks_, options_.borderDensityInverse);
    }
    fill(field_.sizeMax);
    // Without adaption and without borders every size is the target edge length.
    if (options_.densityRatio != 1.0f || border_) {
        progress_.check();
        if (adapt()) {
            fit();
        }
    }
    applyDensity();
    field_.estimatedQuads = SizeFieldBuilder::estimateQuads(mesh_, field_.sizeU, field_.sizeV);
    if (border_) {
        field_.border = border_->release();
    }
    part_.sizeField = std::move(field_);
}

void Computation::setLimits()
{
    const float h = field_.edgeLength;
    if (exactCount()) {
        // Room to both sides, so that the counting can move the sizes.
        const float w = std::sqrt(options_.densityRatio);
        field_.sizeMin = h / w;
        field_.sizeMax = h * w;
    } else {
        field_.sizeMax = h;
        field_.sizeMin = h / options_.densityRatio;
    }
}

void Computation::fill(float size)
{
    field_.sizeU.assign(faceCount_, size);
    field_.sizeV.assign(faceCount_, size);
}

bool Computation::adapt()
{
    const bool adaptive = options_.densityRatio != 1.0f;
    int round = 0;
    float floor = floorShare * field_.sizeMin;
    float start = field_.sizeMax;
    float pointsPerCircle = options_.pointsPerCircle;
    while (true) {
        if (adaptive) {
            followCurvature(pointsPerCircle, floor);
            if (cuts_.mirrorHalfCount != 0) {
                mirror();
            }
        }
        if (border_) {
            border_->apply(field_.sizeU, field_.sizeV);
        }
        if (!adaptive) {
            return false;
        }
        progress_.check();
        smooth();
        clamp(floor);
        if (border_) {
            progress_.check();
            border_->apply(field_.sizeU, field_.sizeV);
        }
        if (!exactCount()) {
            return true;
        }
        progress_.check();
        // The budget as it was given, not the one raised for the target edge length.
        const float q = static_cast<float>(SizeFieldBuilder::estimateQuads(mesh_, field_.sizeU, field_.sizeV)) /
                        static_cast<float>(part_.quadBudget);
        if (!(q > bandTop)) {
            if (!(wide(q) < bandBottom) || round > 3) {
                return true;
            }
            // Too few quads: smaller sizes, a lower floor and from the second round on more points per turn.
            const float g = (q + 1.0f) * 0.5f;
            start = start * g;
            floor = floor * g;
            const float t = round == 0 ? 1.0f : (round == 1 ? static_cast<float>((wide(q) + 3.0) * 0.25) : g);
            pointsPerCircle = pointsPerCircle / t;
            ++round;
        } else if (round == 0) {
            // Too many quads: larger sizes first, then fewer points per turn.
            start = start * minss((q + 1.0f) * 0.5f, 1.25f);
            round = 1;
        } else if (round == 1) {
            const float g = minss((q + 1.0f) * 0.5f, 1.15f);
            pointsPerCircle = pointsPerCircle / g;
            start = start * g;
            round = 2;
        } else if (round == 2 || round == 3) {
            pointsPerCircle = maxss(5.0f, pointsPerCircle / ((q + q + 1.0f) / 3.0f));
            ++round;
        } else {
            return true;
        }
        fill(start);
    }
}

void Computation::followCurvature(float pointsPerCircle, float floor)
{
    const float chord = static_cast<float>(2.0 * std::sin(halfTurn / wide(pointsPerCircle)));
    // The chord of five segments per turn over the curvature, but not above the floor.
    auto least = [floor](float curvature) { return minss(fiveSegments / curvature, floor); };
    for (std::size_t t = 0; t < faceCount_; ++t) {
        if (inMirrorHalf(t)) {
            continue;
        }
        // The curvature along the stored direction limits the size along it, the curvature across it
        // the size across it; which combed direction the stored one is follows from the quarter
        // turns of the triangle.
        const bool even = (cuts_.rotation[t] & 1) == 0;
        float& across = even ? field_.sizeV[t] : field_.sizeU[t];
        float& along = even ? field_.sizeU[t] : field_.sizeV[t];
        const float kAcross = std::fabs(curvature_.lambda2[t]);
        const float kAlong = std::fabs(curvature_.lambda1[t]);
        if (!(curvature_.weight[t] > trustLimit)) {
            // No reliable direction: one curvature for both sizes.
            if (flatLimit > wide(kAcross)) {
                continue;
            }
            const float k = (kAcross + kAcross + kAlong) / 3.0f;
            const float size = chord / k;
            if (size < across) {
                across = size;
            }
            if (size < along) {
                along = size;
            }
            const float low = least(k);
            if (across < low) {
                across = low;
            }
            if (along < low) {
                along = low;
            }
        } else if (!(flatLimit > wide(kAcross) && flatLimit > wide(kAlong))) {
            const float sizeAcross = chord / kAcross;
            if (sizeAcross < across) {
                across = sizeAcross;
            }
            const float sizeAlong = chord / kAlong;
            if (sizeAlong < along) {
                along = sizeAlong;
            }
            const float lowAlong = least(kAlong);
            const float lowAcross = least(kAcross);
            if (across < lowAcross) {
                across = lowAcross;
            }
            if (along < lowAlong) {
                along = lowAlong;
            }
        }
    }
}

void Computation::mirror()
{
    prep::SymmetryCache* symmetric = part_.mirror.symmetric(mesh_, true);
    if (symmetric == nullptr) {
        return;
    }
    const std::vector<int>& faceMirror = symmetric->faceMirror(mesh_);
    for (std::size_t t = 0; t < faceCount_; ++t) {
        if (!inMirrorHalf(t)) {
            continue;
        }
        const int m = faceMirror[t];
        if (m != -1) {
            field_.sizeU[t] = field_.sizeU[static_cast<std::size_t>(m)];
            field_.sizeV[t] = field_.sizeV[static_cast<std::size_t>(m)];
        }
    }
}

void Computation::smooth()
{
    std::vector<float>& u = field_.sizeU;
    std::vector<float>& v = field_.sizeV;
    smoothing_.mean(u, v);
    smoothing_.minimum(2, u, v);
    smoothing_.mean(u, v);
    progress_.check();
    smoothing_.minimum(3, u, v);
    smoothing_.mean(u, v);
    smoothing_.minimum(2, u, v);
    smoothing_.mean(u, v);
}

void Computation::clamp(float floor)
{
    const float ratio = options_.maxQuadRatio;
    for (std::size_t t = 0; t < faceCount_; ++t) {
        if (inMirrorHalf(t)) {
            continue;
        }
        const bool vSmaller = field_.sizeV[t] < field_.sizeU[t];
        float& low = vSmaller ? field_.sizeV[t] : field_.sizeU[t];
        float& high = vSmaller ? field_.sizeU[t] : field_.sizeV[t];
        if (low < floor) {
            low = floor;
        }
        if (low * ratio < high) {
            high = low * ratio;
        }
        if (high < floor) {
            high = floor;
        }
    }
}

Computation::Unknowns Computation::numberUnknowns(int& count) const
{
    Unknowns unknowns(faceCount_, {-1, -1});
    const SizeField::BorderData* border = border_ ? &border_->data() : nullptr;
    count = 0;
    for (std::size_t t = 0; t < faceCount_; ++t) {
        if (inMirrorHalf(t)) {
            continue;
        }
        if (border == nullptr || !border->fixedU.test(t)) {
            unknowns[t][0] = count++;
        }
        if (border == nullptr || !border->fixedV.test(t)) {
            unknowns[t][1] = count++;
        }
    }
    return unknowns;
}

Computation::Weights Computation::dataWeights() const
{
    Weights weights(faceCount_, {1.0, 1.0});
    const float span = field_.sizeMax - field_.sizeMin;
    // A size keeps the more weight the further curvature brought it down from the largest size.
    auto fromSize = [&](float size) {
        float x = (field_.sizeMax - size) / span;
        x = x < 0.0f ? 0.0f : (1.0f < x ? 1.0f : x);
        return wide((x * dataGain + 1.0f) * dataWeight);
    };
    for (std::size_t t = 0; t < faceCount_; ++t) {
        if (inMirrorHalf(t)) {
            continue;
        }
        if ((curvature_.flags[t] & curvature_flags::reliable) != 0) {
            weights[t] = {fromSize(field_.sizeU[t]), fromSize(field_.sizeV[t])};
        }
        if (!border_ || !border_->inSeam(t)) {
            continue;
        }
        const int c = border_->data().borderCorner[t];
        if (part_.lockedBorderEdges.testBounded(static_cast<std::size_t>(adjacency_.cornerEdge(c)))) {
            continue;
        }
        // Beside a frozen border the size across the border holds.
        const bool alongU = cuts_.curveFollowsU.testBounded(static_cast<std::size_t>(c));
        weights[t][alongU == border_->swapped(t) ? 0 : 1] *= borderWeight;
    }
    const std::vector<float>& doubleArea = mesh_.triangleGeometry().doubleAreas();
    const float area = mesh_.triangleGeometry().totalArea();
    const float perArea = static_cast<float>(mesh_.faceCount()) / (area + area);
    for (std::size_t t = 0; t < faceCount_; ++t) {
        if (inMirrorHalf(t)) {
            continue;
        }
        const double g = wide(doubleArea[t] * perArea);
        weights[t][0] *= g;
        weights[t][1] *= g;
    }
    return weights;
}

void Computation::addEdgeRow(const Unknowns& unknowns, const Weights& weights, std::size_t t0, int p, std::size_t t1,
                             int q, Rows& rows) const
{
    // Size p of t0 equals size q of t1. A fixed size goes to the target; a
    // row between two fixed sizes says nothing and is left out with its target.
    solver::ListMatrix::List entries;
    double target = 0.0;
    auto term = [&](std::size_t t, int component, double sign) {
        const int unknown = unknowns[t][static_cast<std::size_t>(component)];
        if (unknown != -1) {
            entries.push_back({unknown, sign});
        } else {
            target = target - sign * wide(size(t, component));
        }
    };
    term(t0, p, 1.0);
    term(t1, q, -1.0);
    if (entries.empty()) {
        return;
    }
    // The first corner of an inner edge lies in the triangle with the smaller number.
    assert(entries.size() < 2 || entries[0].index < entries[1].index);
    const double w0 = weights[t0][static_cast<std::size_t>(p)];
    const double w1 = weights[t1][static_cast<std::size_t>(q)];
    rows.entries.push_back(std::move(entries));
    rows.targets.push_back(target);
    rows.weights.push_back(maxsd(edgeWeightMin, ((w0 + w1) * 0.5) * edgeWeightFactor));
}

void Computation::fit()
{
    int unknownCount = 0;
    const Unknowns unknowns = numberUnknowns(unknownCount);
    const Weights weights = dataWeights();
    progress_.check();

    Rows rows;
    // Every free size stays near its value.
    for (std::size_t t = 0; t < faceCount_; ++t) {
        if (inMirrorHalf(t)) {
            continue;
        }
        for (int component = 0; component < 2; ++component) {
            const int unknown = unknowns[t][static_cast<std::size_t>(component)];
            if (unknown != -1) {
                rows.entries.push_back({{unknown, 1.0}});
                rows.targets.push_back(wide(size(t, component)));
                rows.weights.push_back(weights[t][static_cast<std::size_t>(component)]);
            }
        }
    }
    progress_.check();
    // Across every inner edge the sizes along the same direction agree; across a
    // curve edge only the sizes along the curve.
    for (int e = 0; e < adjacency_.edgeCount(); ++e) {
        const auto faces = adjacency_.innerEdgeFaces(e);
        if (!faces) {
            continue;
        }
        const auto t0 = static_cast<std::size_t>((*faces)[0]);
        const auto t1 = static_cast<std::size_t>((*faces)[1]);
        if (inMirrorHalf(t0) || inMirrorHalf(t1)) {
            continue;
        }
        const auto edge = static_cast<std::size_t>(e);
        const bool even = (cuts_.residualJump[edge] & 1) == 0;
        const std::array<std::array<int, 2>, 2> pairs = even ? std::array<std::array<int, 2>, 2>{{{0, 0}, {1, 1}}}
                                                             : std::array<std::array<int, 2>, 2>{{{0, 1}, {1, 0}}};
        if (cuts_.curveEdgeKind[edge] != curve_edge::none) {
            const int first = adjacency_.edgeCorners(e)[0];
            const auto& pair = pairs[cuts_.curveFollowsU.testBounded(static_cast<std::size_t>(first)) ? 0 : 1];
            addEdgeRow(unknowns, weights, t0, pair[0], t1, pair[1], rows);
        } else {
            for (const auto& pair : pairs) {
                addEdgeRow(unknowns, weights, t0, pair[0], t1, pair[1], rows);
            }
        }
    }
    progress_.check();

    solver::ListMatrix matrix(static_cast<int>(rows.entries.size()), unknownCount);
    for (std::size_t r = 0; r < rows.entries.size(); ++r) {
        matrix.list(static_cast<int>(r)) = std::move(rows.entries[r]);
    }
    solver::NormalEquations equations;
    if (equations.assemble(matrix, rows.weights, rows.targets, solver::NormalEquations::Scope::All, progress_) !=
        solver::NormalEquations::Result::Ok) {
        return;
    }
    progress_.check();
    solver::DirectSolver solver;
    std::vector<double> x;
    // Without a solution the field stays as smoothing and clamping left it.
    if (!solver.solve(equations.matrix(), false, equations.rightHandSide(), progress_, x)) {
        return;
    }
    for (std::size_t t = 0; t < faceCount_; ++t) {
        if (unknowns[t][0] != -1) {
            field_.sizeU[t] = static_cast<float>(x[static_cast<std::size_t>(unknowns[t][0])]);
        }
        if (unknowns[t][1] != -1) {
            field_.sizeV[t] = static_cast<float>(x[static_cast<std::size_t>(unknowns[t][1])]);
        }
    }
    if (cuts_.mirrorHalfCount != 0) {
        mirror();
    }
}

void Computation::applyDensity()
{
    std::optional<std::vector<float>> factors;
    if (options_.valueDensity) {
        factors = DensityFactors::fromValues(mesh_);
    }
    if (!factors && options_.colorDensity) {
        factors = DensityFactors::fromColors(mesh_);
    }
    if (factors) {
        DensityFactors::apply(mesh_, *factors, field_.sizeU, field_.sizeV);
    }
}

} // namespace

void SizeFieldBuilder::compute(parts::Part& part, const Options& options, const run::Progress& progress)
{
    part.sizeField.reset();
    if (!part.curvature || !part.crossField || !part.cutGraph) {
        return;
    }
    Computation(part, options, progress).run();
}

int SizeFieldBuilder::estimateQuads(const Mesh& mesh, std::span<const float> sizeU, std::span<const float> sizeV)
{
    const TriangleGeometry& geometry = mesh.triangleGeometry();
    const std::vector<float>& doubleArea = geometry.doubleAreas();
    const float area = geometry.totalArea();
    float quadArea = 0.0f;
    float quads = 0.0f;
    for (std::size_t t = 0; t < sizeU.size(); ++t) {
        quadArea += (doubleArea[t] * sizeU[t]) * sizeV[t];
        quads += doubleArea[t] / (sizeV[t] * sizeU[t]);
    }
    quadArea = quadArea * 0.5f;
    quads = quads * 0.5f;
    return truncToInt((area / (quadArea / area) + quads) * 0.5f);
}

} // namespace areko::remesher::field
