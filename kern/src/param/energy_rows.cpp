#include "energy_rows.h"

#include "mesh/vec3.h"
#include "param/constraint_terms.h"

#include <array>

namespace areko::remesher::param {

using namespace areko::remesher::mesh;

namespace {

// Weight factor of the rows that ask for a growth of 1.
constexpr float unitRowWeight = 0.55f;
// Stands in for the cross product of a triangle without area.
constexpr float flatTriangle = 1e-5f;

// Coefficients of one row at the three corners: the gradient of the hat
// function of each corner of the flat triangle (0, 0), (x1, 0), (x2, y2)
// along `direction`, times twice the triangle area and `scale`.
[[nodiscard]] std::array<double, 3> coefficients(const std::array<float, 2>& direction, double x1, double x2, double y2,
                                                 double scale)
{
    const double dx = x2 - x1;
    const double first = wide(direction[0]);
    const double second = wide(direction[1]);
    return {(dx * second - first * y2) * scale, (first * y2 - x2 * second) * scale, second * (x1 * scale)};
}

// Terms of a row at the variables `variables`, corners in order; a
// coefficient that is exactly 0 gets no term.
[[nodiscard]] TermList rowOf(const std::array<int, 3>& variables, const std::array<double, 3>& values)
{
    TermList row;
    for (std::size_t i = 0; i < 3; ++i) {
        if (values[i] != 0.0) {
            row.push_back({variables[i], values[i]});
        }
    }
    return row;
}

} // namespace

EnergyRows::Result EnergyRows::build(const SystemContext& context)
{
    const parts::Part& part = context.part();
    const Topology& topology = context.mesh().topology();
    const std::vector<Vec3f>& points = context.mesh().points();
    const field::CrossField& frames = *part.crossField;
    const field::CutGraph& graph = context.graph();
    const std::size_t faces = static_cast<std::size_t>(context.mesh().faceCount());
    const BitArray& mirrorHalf = graph.mirrorHalf;

    const std::size_t rowCount = 4 * faces - (context.symmetry() ? 4 * static_cast<std::size_t>(graph.mirrorHalfCount) : 0);
    Result out;
    out.rows.assign(rowCount, TermList{});
    out.targets.assign(rowCount, 0.0);
    out.weights.assign(rowCount, 1.0);
    out.triangleRowStart.reserve(faces + 1);
    out.triangleRowStart.push_back(0);

    std::size_t z = 0;
    for (std::size_t f = 0; f < faces; ++f) {
        if (context.symmetry() && mirrorHalf.testBounded(f)) {
            out.triangleRowStart.push_back(static_cast<int>(z));
            continue;
        }
        const int base = 3 * static_cast<int>(f);
        const std::uint8_t start = frames.startCorner[f];
        const std::array<int, 3> corners = start == 0   ? std::array<int, 3>{base, base + 1, base + 2}
                                           : start == 1 ? std::array<int, 3>{base + 1, base + 2, base}
                                                        : std::array<int, 3>{base + 2, base, base + 1};
        const Vec3f a = points[static_cast<std::size_t>(topology.vertex(corners[0]))];
        const Vec3f e1 = points[static_cast<std::size_t>(topology.vertex(corners[1]))] - a;
        const Vec3f e2 = points[static_cast<std::size_t>(topology.vertex(corners[2]))] - a;
        float doubleArea = length(cross(e1, e2));
        if (doubleArea == 0.0f) {
            doubleArea = flatTriangle;
        }
        const double area = wide(doubleArea);
        const double sizeU = wide(context.sizeAlong(f, 0));
        const double sizeV = wide(context.sizeAlong(f, 1));
        const double x1 = wide(frames.referenceLength[f]);
        const double x2 = wide(dot(frames.frameE1[f], e2));
        const double y2 = wide(dot(e2, frames.frameE2[f]));
        const double scaleU = sizeU / area;
        const double scaleV = sizeV / area;
        const std::array<float, 2>& alongFirst = graph.combedU[f];
        const std::array<float, 2>& alongSecond = graph.combedV[f];

        std::array<int, 3> u;
        std::array<int, 3> v;
        for (std::size_t i = 0; i < 3; ++i) {
            u[i] = 2 * context.copyPoint(corners[i]);
            v[i] = u[i] + 1;
        }
        const double unitWeight = wide(doubleArea * unitRowWeight);
        out.rows[z] = rowOf(u, coefficients(alongFirst, x1, x2, y2, scaleU));
        out.targets[z] = 1.0;
        out.weights[z] = unitWeight;
        out.rows[z + 1] = rowOf(u, coefficients(alongSecond, x1, x2, y2, scaleU));
        out.weights[z + 1] = area;
        out.rows[z + 2] = rowOf(v, coefficients(alongSecond, x1, x2, y2, scaleV));
        out.targets[z + 2] = 1.0;
        out.weights[z + 2] = unitWeight;
        out.rows[z + 3] = rowOf(v, coefficients(alongFirst, x1, x2, y2, scaleV));
        out.weights[z + 3] = area;
        z += 4;
        out.triangleRowStart.push_back(static_cast<int>(z));
    }
    for (TermList& row : out.rows) {
        if (!row.empty()) {
            ConstraintTerms::normalize(row);
        }
    }
    return out;
}

} // namespace areko::remesher::param
