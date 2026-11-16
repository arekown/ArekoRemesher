#include "region_bridges.h"

#include "util/introsort.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>

namespace areko::remesher::param {

namespace {

// Cost of an edge of a chain between two sources: below every real cost.
constexpr float meetingChains = -1e-7f;

} // namespace

void RegionBridges::add(const RegionGrowth& growth, const UvOutput::Coordinates& uv, const mesh::Adjacency& copy)
{
    const std::size_t firstOfCoordinate = list_.size();
    std::vector<float> costs;
    for (int e = 0; e < copy.edgeCount(); ++e) {
        const std::array<int, 2> ends = copy.edgeVertices(e);
        const int ra = growth.region(ends[0]);
        const int rb = growth.region(ends[1]);
        if (ra == -1 || rb == -1 || ra == rb) {
            continue;
        }
        float cost = (growth.distance(ends[0]) + growth.distance(ends[1])) + gridDistance(uv, ends[0], ends[1]);
        int edge = e;
        if (growth.onChain(e) && growth.distance(ends[0]) == 0.0f && growth.distance(ends[1]) == 0.0f) {
            cost = meetingChains;
            edge = -1;
        }
        const int low = std::min(ra, rb);
        const int high = std::max(ra, rb);
        std::size_t found = firstOfCoordinate;
        while (found < list_.size() && !(list_[found].low == low && list_[found].high == high)) {
            ++found;
        }
        if (found == list_.size()) {
            Bridge bridge;
            bridge.low = low;
            bridge.high = high;
            bridge.lowIsChain = low < growth.chainRegions();
            bridge.highIsChain = high < growth.chainRegions();
            bridge.edge = edge;
            bridge.coordinate = growth.coordinate();
            list_.push_back(bridge);
            costs.push_back(cost);
        } else if (cost < costs[found - firstOfCoordinate]) {
            list_[found].edge = edge;
            costs[found - firstOfCoordinate] = cost;
        }
    }

    const auto c = static_cast<std::size_t>(growth.coordinate());
    for (std::size_t i = firstOfCoordinate; i < list_.size(); ++i) {
        Bridge& bridge = list_[i];
        if (bridge.edge == -1) {
            continue;
        }
        const std::array<int, 2> ends = copy.edgeVertices(bridge.edge);
        const int rootA = growth.root(ends[0]);
        const int rootB = growth.root(ends[1]);
        bridge.first = growth.coordinate() + 2 * rootA;
        bridge.second = growth.coordinate() + 2 * rootB;
        bridge.span = std::fabs(uv[static_cast<std::size_t>(rootB)][c] - uv[static_cast<std::size_t>(rootA)][c]);
    }
}

std::vector<int> RegionBridges::bySpan() const
{
    std::vector<int> order(list_.size());
    std::iota(order.begin(), order.end(), 0);
    util::Introsort::sort(order.begin(), order.end(), [this](int a, int b) {
        return list_[static_cast<std::size_t>(a)].span < list_[static_cast<std::size_t>(b)].span;
    });
    return order;
}

} // namespace areko::remesher::param
