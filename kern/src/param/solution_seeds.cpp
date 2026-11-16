#include "solution_seeds.h"

#include "mesh/bit_array.h"
#include "mesh/vec3.h"

#include <utility>

namespace areko::remesher::param {

using namespace areko::remesher::mesh;

void SolutionSeeds::addBorderPoints(const SystemContext& context, bool freezeBorders, Solution& solution)
{
    if (!solution.singularities.points.empty() || !(freezeBorders || solution.seeds.empty())) {
        return;
    }
    appendBorderPoints(context.copy(), context.mesh().pointCount(), solution.uv, solution.seeds);
}

void SolutionSeeds::appendBorderPoints(const field::CutCopy& copy, int pointCount,
                                       std::span<const std::array<float, 2>> uv, std::vector<int>& seeds)
{
    BitArray listed(static_cast<std::size_t>(pointCount));
    for (const int point : seeds) {
        listed.setBounded(static_cast<std::size_t>(point));
    }
    const std::vector<int>& origin = copy.pointOrigins();
    for (const int b : copy.adjacency().boundaryPoints()) {
        const std::array<float, 2>& at = uv[static_cast<std::size_t>(b)];
        if (!isWholeNumber(at[0]) || !isWholeNumber(at[1])) {
            continue;
        }
        const auto point = static_cast<std::size_t>(origin[static_cast<std::size_t>(b)]);
        if (!listed.test(point)) {
            listed.set(point);
            seeds.push_back(static_cast<int>(point));
        }
    }
}

void SolutionSeeds::mirror(const SystemContext& context, const prep::SymmetryCache& pairs, Solution& solution)
{
    const field::CutGraph& graph = context.graph();
    if (graph.mirrorHalfCount == 0) {
        return;
    }
    // Points that touch only triangles of the mirrored half.
    const Topology& topology = context.mesh().topology();
    BitArray mirroredOnly(static_cast<std::size_t>(context.mesh().pointCount()), true);
    for (std::size_t t = 0; t < graph.mirrorHalf.size(); ++t) {
        if (graph.mirrorHalf.test(t)) {
            continue;
        }
        for (int k = 0; k < 3; ++k) {
            mirroredOnly.setBounded(static_cast<std::size_t>(topology.vertex(3 * static_cast<int>(t) + k)), false);
        }
    }
    const std::vector<int>& pointMirror = pairs.pointMirror();
    // The mirror point of `point` when it is another point, else -1.
    auto partner = [&](int point) {
        const auto p = static_cast<std::size_t>(point);
        const int m = p < pointMirror.size() ? pointMirror[p] : -1;
        return m != point ? m : -1;
    };

    field::Singularities kept;
    const field::Singularities& all = solution.singularities;
    auto keep = [&](int point, std::size_t j) {
        kept.points.push_back(point);
        kept.index.push_back(all.index[j]);
        kept.fullTurn.push_back(all.fullTurn[j]);
        kept.indexSum += all.index[j];
    };
    for (std::size_t j = 0; j < all.points.size(); ++j) {
        const int s = all.points[j];
        if (mirroredOnly.testBounded(static_cast<std::size_t>(s))) {
            continue;
        }
        keep(s, j);
        if (const int m = partner(s); m != -1) {
            keep(m, j);
        }
    }
    std::vector<int> seeds;
    for (const int q : solution.seeds) {
        if (mirroredOnly.testBounded(static_cast<std::size_t>(q))) {
            continue;
        }
        seeds.push_back(q);
        if (const int m = partner(q); m != -1) {
            seeds.push_back(m);
        }
    }
    solution.singularities = std::move(kept);
    solution.seeds = std::move(seeds);
}

} // namespace areko::remesher::param
