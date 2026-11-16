#include "bridge_rows.h"

#include "mesh/vec3.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace areko::remesher::param {

using namespace areko::remesher::mesh;

namespace {

// What a bridge asks for.
enum class Demand {
    Soft,          // the present difference, weighted
    SoftAndCount,  // the same, and the difference is the count of the bridge
    Unit,          // the difference is exactly 1
};

struct Rows {
    Demand demand = Demand::Soft;
    float weight = 1.0f;
    double target = 0.0;
};

// Below this span two sources are a single grid step apart.
constexpr float singleStep = 1.35f;

[[nodiscard]] double unitTowards(double target)
{
    return 0.0 > target ? -1.0 : 1.0;
}

// `along` and `across` are the distances of the two sources in the
// coordinate of the bridge and in the other one, `target` the signed
// difference along. A distance that is not a number leads to a soft row.
[[nodiscard]] Rows demandOf(const Bridge& bridge, float along, float across, double target)
{
    if (bridge.lowIsChain && bridge.highIsChain) {
        if (!(singleStep > along)) {
            return {5.0f > along && along > across ? Demand::SoftAndCount : Demand::Soft, 1.0f, target};
        }
        if (!(along + along > across)) {
            return {Demand::Soft, 1.0f, target};
        }
        return {Demand::Unit, 1.0f, unitTowards(target)};
    }
    if (bridge.lowIsChain || bridge.highIsChain) {
        if (!(along + along > across)) {
            return {Demand::Soft, 1.0f, target};
        }
        if (!(singleStep > along)) {
            return {4.5f > along ? Demand::SoftAndCount : Demand::Soft, 1.0f, target};
        }
        return {Demand::SoftAndCount, 1.0f, unitTowards(target)};
    }
    // Two singular points: the weight falls from 100 for a short straight
    // bridge to 1 for a long one, and to 0.01 for one that runs mostly
    // across the coordinate.
    float weight = 0.01f;
    if (0.2 * wide(across) > wide(along)) {
        return {Demand::Soft, weight, target};
    }
    if (across > along) {
        weight = ((across / along - 0.2f) * 0.09f) / 0.8f + 0.01f;
    } else if (1.0f > along) {
        weight = 100.0f;
    } else if (2.5f > along) {
        weight = ((along - 1.0f) * -90.0f) / 1.5f + 100.0f;
    } else if (5.0f > along) {
        weight = ((along - 2.5f) * -9.0f) / 2.5f + 10.0f;
    } else {
        weight = 1.0f;
    }
    return {Demand::Soft, weight, target};
}

// How many grid steps of `coordinate` the combed field counts along a path
// of copy points: per edge the step in space along the combed direction over
// the size, averaged over the triangles of the edge (float throughout).
// Positive when the coordinate grows along the path.
[[nodiscard]] float fieldAlong(const SystemContext& context, std::span<const int> path, int coordinate)
{
    const parts::Part& part = context.part();
    const Adjacency& copy = context.copy().adjacency();
    const std::vector<int>& origin = context.copy().pointOrigins();
    const std::vector<Vec3f>& points = context.mesh().points();
    const field::CrossField& frames = *part.crossField;
    const std::vector<std::array<float, 2>>& combed = coordinate == 0 ? context.graph().combedU
                                                                       : context.graph().combedV;
    float sum = 0.0f;
    for (std::size_t i = 0; i + 1 < path.size(); ++i) {
        const int e = copy.edgeBetween(path[i], path[i + 1]);
        if (e == -1) {
            continue;
        }
        const Vec3f step = points[static_cast<std::size_t>(origin[static_cast<std::size_t>(path[i + 1])])]
                           - points[static_cast<std::size_t>(origin[static_cast<std::size_t>(path[i])])];
        const std::span<const int> corners = copy.edgeCorners(e);
        float steps = 0.0f;
        for (const int corner : corners) {
            const auto t = static_cast<std::size_t>(corner / 3);
            steps = steps + dot(frames.inSpace(t, combed[t]), step) / context.sizeAlong(t, coordinate);
        }
        sum = steps / static_cast<float>(corners.size()) + sum;
    }
    return sum;
}

// The way from the source of the first end of the bridge edge over the edge
// to the source of its second end.
[[nodiscard]] std::vector<int> wayBetweenSources(const Bridge& bridge, const RegionGrowth& growth,
                                                 const Adjacency& copy)
{
    const std::array<int, 2> ends = copy.edgeVertices(bridge.edge);
    std::vector<int> head = growth.pathToRoot(ends[0]);
    std::vector<int> tail = growth.pathToRoot(ends[1]);
    // The first variable of a bridge belongs to the source of its first end
    // (RegionBridges), so the two paths are never exchanged here.
    if (bridge.first / 2 != head.back()) {
        head.swap(tail);
    }
    std::vector<int> way(head.rbegin(), head.rend());
    way.insert(way.end(), tail.begin(), tail.end());
    return way;
}

} // namespace

void BridgeRows::add(const SystemContext& context, const ConstraintSet& constraints, const UvOutput::Coordinates& uv,
                     std::span<const RegionGrowth, 2> growth, const RegionBridges& bridges, SoftSystem& system)
{
    const Adjacency& copy = context.copy().adjacency();
    for (const int number : bridges.bySpan()) {
        const Bridge& bridge = bridges.list()[static_cast<std::size_t>(number)];
        if (bridge.first == -1) {
            continue;
        }
        const int a = bridge.first / 2;
        const int b = bridge.second / 2;
        if (context.symmetric() && (context.mirrored(a) || context.mirrored(b))) {
            continue;
        }
        const auto c = static_cast<std::size_t>(bridge.coordinate);
        const std::array<float, 2>& from = uv[static_cast<std::size_t>(a)];
        const std::array<float, 2>& to = uv[static_cast<std::size_t>(b)];
        const float along = to[c] - from[c];
        const float across = std::fabs(to[1 - c] - from[1 - c]);
        Rows rows = demandOf(bridge, std::fabs(along), across, wide(along));
        if (rows.demand != Demand::Soft) {
            // The field decides on which side the difference is positive.
            const std::vector<int> way = wayBetweenSources(bridge, growth[c], copy);
            if (0.0f > along * fieldAlong(context, way, bridge.coordinate)) {
                rows.target = -rows.target;
            }
        }

        if (rows.demand == Demand::Unit) {
            TermList row{{bridge.second, 1.0}, {bridge.first, -1.0}};
            solver::ListMatrix::order(row, solver::ListMatrix::SameIndex::Separate);
            constraints.removeAnchorTerms(row);
            if (!row.empty()) {
                system.addHard(std::move(row), rows.target);
            }
            continue;
        }
        TermList soft{{bridge.second, 1.0}, {bridge.first, -1.0}};
        constraints.removeAnchorTerms(soft);
        system.addSoft(std::move(soft), rows.target, wide(rows.weight));
        if (rows.demand == Demand::SoftAndCount) {
            const double sign = rows.target > 0.0 ? 1.0 : -1.0;
            TermList count{{bridge.second, sign}, {bridge.first, -sign}};
            constraints.removeAnchorTerms(count);
            solver::ListMatrix::order(count, solver::ListMatrix::SameIndex::Separate);
            count.push_back({system.bridgeColumn(number), -1.0});
            system.addHard(std::move(count), 0.0);
        }
    }
}

} // namespace areko::remesher::param
