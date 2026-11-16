#include "uv_output.h"

#include "mesh/vec3.h"

namespace areko::remesher::param {

using namespace areko::remesher::mesh;

std::vector<float> UvOutput::variables(const System& system, std::span<const double> x)
{
    std::vector<float> values(static_cast<std::size_t>(system.variableCount), 0.0f);
    for (std::size_t i = 0; i < values.size(); ++i) {
        const int unknown = system.variableToUnknown[i];
        if (unknown >= 0) {
            values[i] = static_cast<float>(x[static_cast<std::size_t>(unknown)]);
        }
    }
    const std::vector<Constraint>& conditions = system.constraints.list();
    for (std::size_t b = conditions.size(); b-- > 0;) {
        const Constraint& condition = conditions[b];
        if (condition.pivot < 0) {
            continue;
        }
        float value = 0.0f;
        if (!condition.terms.empty()) {
            double sum = 0.0;
            for (const Term& term : condition.terms) {
                sum = sum + wide(values[static_cast<std::size_t>(term.index)]) * term.value;
            }
            value = static_cast<float>(sum);
        }
        values[static_cast<std::size_t>(condition.pivot)] = value;
    }
    return values;
}

UvOutput::Coordinates UvOutput::coordinates(const System& system, std::span<const double> x)
{
    const std::vector<float> values = variables(system, x);
    Coordinates uv(static_cast<std::size_t>(system.copyPointCount));
    for (std::size_t w = 0; w < uv.size(); ++w) {
        uv[w] = {values[2 * w], values[2 * w + 1]};
    }
    return uv;
}

std::array<int, 3> UvOutput::mirroredCorners(const Topology& topology, std::span<const int> pointMirror,
                                             std::span<const int> faceMirror, int t)
{
    static constexpr std::array<int, 3> none{-1, -1, -1};
    const auto face = static_cast<std::size_t>(t);
    if (face >= faceMirror.size() || faceMirror[face] < 0) {
        return none;
    }
    const int first = 3 * faceMirror[face];
    std::array<int, 3> mirror{};
    for (int j = 0; j < 3; ++j) {
        const auto point = static_cast<std::size_t>(topology.vertex(3 * t + j));
        mirror[static_cast<std::size_t>(j)] = point < pointMirror.size() ? pointMirror[point] : -1;
    }
    const int p0 = topology.vertex(first);
    const int p1 = topology.vertex(first + 1);
    const int p2 = topology.vertex(first + 2);
    if (p0 == mirror[0]) {
        return p2 == mirror[1] && p1 == mirror[2] ? std::array<int, 3>{first, first + 2, first + 1} : none;
    }
    if (p0 == mirror[1]) {
        return p1 == mirror[0] && p2 == mirror[2] ? std::array<int, 3>{first + 1, first, first + 2} : none;
    }
    if (p0 == mirror[2]) {
        return p2 == mirror[0] && p1 == mirror[1] ? std::array<int, 3>{first + 2, first + 1, first} : none;
    }
    return none;
}

void UvOutput::mirror(const SystemContext& context, prep::SymmetryCache& pairs, Coordinates& uv)
{
    const Topology& topology = context.mesh().topology();
    const std::vector<int>& faceMirror = pairs.faceMirror(context.mesh());
    const std::vector<int>& pointMirror = pairs.pointMirror();
    context.graph().mirrorHalf.forEachSet([&](std::size_t t) {
        const int base = 3 * static_cast<int>(t);
        const std::array<int, 3> source = mirroredCorners(topology, pointMirror, faceMirror, static_cast<int>(t));
        for (int k = 0; k < 3; ++k) {
            const int corner = source[static_cast<std::size_t>(k)];
            // A corner without mirror corner has nothing to take over.
            if (corner < 0) {
                continue;
            }
            const std::array<float, 2> from = uv[static_cast<std::size_t>(context.copyPoint(corner))];
            uv[static_cast<std::size_t>(context.copyPoint(base + k))] = {-from[0], from[1]};
        }
    });
}

} // namespace areko::remesher::param
