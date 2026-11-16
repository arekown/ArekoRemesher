#include "reweighting.h"

#include "field/smoothing_neighbours.h"
#include "mesh/vec3.h"
#include "param/uv_output.h"

#include <array>
#include <cmath>
#include <vector>

namespace areko::remesher::param {

using namespace areko::remesher::mesh;

namespace {

// Surcharge of a triangle that is turned over, and the cap of every factor.
constexpr double flippedSurcharge = 5.0;
constexpr double largestFactor = 10.0;
// A grid step counts as too long from this multiple of the size on, and its
// surcharge is full at the second multiple.
constexpr float stepStart = 1.8f;
constexpr float stepFull = 4.0f;
// The angle between the grid directions counts from this far off a right
// angle on; the surcharge grows over the span up to 90 degrees off.
constexpr double angleStart = degD(15.0);
constexpr double angleSpan = degD(75.0);

// Gradients of the two grid coordinates over a triangle, as vectors in space
// (float throughout). One over the length of a gradient is the length of
// one grid step in the triangle.
struct Gradients {
    Vec3f u;
    Vec3f v;
};

[[nodiscard]] Gradients gradientsOf(Vec3f normal, const std::array<Vec3f, 3>& p,
                                    const std::array<std::array<float, 2>, 3>& uv)
{
    Vec3f e = p[1] - p[0];
    Vec3f w = cross(normal, e);
    const float length = normalizeWithLength(e);
    normalizeInPlace(w);
    const Vec3f d = p[2] - p[0];
    const float a = dot(d, e);
    const float b = dot(d, w);
    const float q = length * b;
    const float k = a / q - 1.0f / b;
    const float alongU = (uv[1][0] - uv[0][0]) / length;
    const float acrossU = (uv[0][0] * k - (uv[1][0] * a) / q) + uv[2][0] / b;
    const float alongV = (uv[1][1] - uv[0][1]) / length;
    const float acrossV = (k * uv[0][1] - (a * uv[1][1]) / q) + uv[2][1] / b;
    return {{acrossU * w.x + alongU * e.x, acrossU * w.y + alongU * e.y, acrossU * w.z + alongU * e.z},
            {w.x * acrossV + e.x * alongV, w.y * acrossV + e.y * alongV, acrossV * w.z + e.z * alongV}};
}

// Normalises both gradients, returns their lengths and the surcharge for
// their angle (0 when they are close enough to a right angle). The angle is
// taken between the normalised gradients and compared in radians.
[[nodiscard]] double angleSurcharge(Gradients& g, float& lengthU, float& lengthV)
{
    lengthU = normalizeWithLength(g.u);
    lengthV = normalizeWithLength(g.v);
    const double cosine = wide(dot(g.u, g.v));
    float offRight = field::quarterTurnF;
    if (!(cosine <= -1.0) && !(cosine >= 1.0)) {
        offRight = std::fabs(static_cast<float>(std::acos(cosine)) - field::quarterTurnF);
    }
    if (!(wide(offRight) > angleStart)) {
        return 0.0;
    }
    return wide(static_cast<float>((wide(offRight) - angleStart) / angleSpan) * 3.0f + 1.0f);
}

// Adds the surcharge for a grid step of length `step` where `size` is asked for.
void addStepSurcharge(double& factor, float step, float size)
{
    const float from = size * stepStart;
    if (step > from) {
        const float full = size * stepFull;
        const double surcharge = step > full ? 3.0 : wide((((step - from) * 3.0f) / (full - from)) + 0.0f);
        factor = surcharge + factor;
    }
}

// One factor per row, before averaging.
[[nodiscard]] std::vector<double> rowFactors(const SystemContext& context, const System& system,
                                             const UvOutput::Coordinates& uv)
{
    const Topology& topology = context.mesh().topology();
    const std::vector<Vec3f>& points = context.mesh().points();
    const std::vector<Vec3f>& normals = context.mesh().triangleGeometry().normals();
    const std::size_t faces = static_cast<std::size_t>(context.mesh().faceCount());

    std::vector<double> factors(system.weights.size(), 1.0);
    for (std::size_t t = 0; t < faces; ++t) {
        const auto r = static_cast<std::size_t>(system.triangleRowStart[t]);
        if (system.triangleRowStart[t] == system.triangleRowStart[t + 1]) {
            continue;
        }
        const int base = 3 * static_cast<int>(t);
        std::array<Vec3f, 3> p;
        std::array<std::array<float, 2>, 3> c;
        for (int i = 0; i < 3; ++i) {
            p[static_cast<std::size_t>(i)] = points[static_cast<std::size_t>(topology.vertex(base + i))];
            c[static_cast<std::size_t>(i)] = uv[static_cast<std::size_t>(context.copyPoint(base + i))];
        }
        Gradients g = gradientsOf(normals[t], p, c);
        float lengthU = 0.0f;
        float lengthV = 0.0f;
        const double forAngle = angleSurcharge(g, lengthU, lengthV);
        const float stepU = 1.0f / lengthU;
        const float stepV = 1.0f / lengthV;

        const float doubleArea = (c[2][1] - c[0][1]) * (c[1][0] - c[0][0]) - (c[2][0] - c[0][0]) * (c[1][1] - c[0][1]);
        if (doubleArea < 0.0f) {
            for (std::size_t i = r; i < r + 4; ++i) {
                factors[i] += flippedSurcharge;
            }
        }
        addStepSurcharge(factors[r], stepU, context.sizeAlong(t, 0));
        addStepSurcharge(factors[r + 2], stepV, context.sizeAlong(t, 1));
        if (forAngle != 0.0) {
            factors[r + 1] += forAngle;
            factors[r + 3] = forAngle + factors[r + 3];
        }
        for (std::size_t i = r; i < r + 4; ++i) {
            if (factors[i] > largestFactor) {
                factors[i] = largestFactor;
            }
        }
    }
    return factors;
}

// Averages the factors of every triangle with those of its neighbours,
// `passes` times; each pass reads the factors of the pass before.
void smooth(const SystemContext& context, const System& system, std::vector<double>& factors, int passes)
{
    if (passes < 1) {
        return;
    }
    const field::SmoothingNeighbours neighbours(context.adjacency(), context.graph(),
                                                field::SmoothingNeighbours::Limit::Creases);
    const std::vector<int>& rowStart = system.triangleRowStart;
    const int faces = context.mesh().faceCount();
    auto firstRow = [&](int triangle) { return rowStart[static_cast<std::size_t>(triangle)]; };

    std::vector<double> next(factors.size(), 0.0);
    for (int pass = 0; pass < passes; ++pass) {
        for (int t = 0; t < faces; ++t) {
            const int row = firstRow(t);
            if (firstRow(t + 1) != row + 4) {
                continue;
            }
            const auto r = static_cast<std::size_t>(row);
            std::array<double, 4> sum{factors[r], factors[r + 1], factors[r + 2], factors[r + 3]};
            float countU = 1.0f;
            float countV = 1.0f;
            for (int k = 0; k < 3; ++k) {
                const auto across = neighbours.across(3 * t + k);
                // Only a neighbour with rows counts.
                if (!across || firstRow(across->triangle) == firstRow(across->triangle + 1)) {
                    continue;
                }
                // The rows the neighbour gives to the first and to the second direction.
                const auto s = static_cast<std::size_t>(firstRow(across->triangle));
                const std::size_t forU = across->odd ? s + 2 : s;
                const std::size_t forV = across->odd ? s : s + 2;
                if (across->take.u) {
                    countU += 1.0f;
                    sum[0] += factors[forU];
                    sum[1] += factors[forU + 1];
                }
                if (across->take.v) {
                    countV += 1.0f;
                    sum[2] += factors[forV];
                    sum[3] += factors[forV + 1];
                }
            }
            next[r] = sum[0] / wide(countU);
            next[r + 1] = sum[1] / wide(countU);
            next[r + 2] = sum[2] / wide(countV);
            next[r + 3] = sum[3] / wide(countV);
        }
        // Rows of no triangle are never written: they carry what the other
        // buffer held, their factor after an even number of passes and 0
        // after an odd one.
        factors.swap(next);
    }
}

} // namespace

void Reweighting::apply(const SystemContext& context, System& system, std::span<const double> x,
                        const run::Progress& progress)
{
    std::vector<double> factors = rowFactors(context, system, UvOutput::coordinates(system, x));
    smooth(context, system, factors, smoothingPasses);
    for (std::size_t i = 0; i < system.weights.size(); ++i) {
        system.weights[i] = system.weights[i] * factors[i];
    }
    (void)system.normal.assemble(system.rows, system.weights, system.targets,
                                 solver::NormalEquations::Scope::KeepTranspose, progress);
}

} // namespace areko::remesher::param
