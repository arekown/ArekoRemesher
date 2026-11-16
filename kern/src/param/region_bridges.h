#pragma once

#include "param/region_growth.h"

#include <vector>

namespace areko::remesher::param {

// Two neighbouring regions of one coordinate and the cheapest edge between
// them: the sources of the two regions lie `span` grid lines apart.
struct Bridge {
    // The two regions, low < high.
    int low = -1;
    int high = -1;
    // Whether the regions belong to curve chains (otherwise to singular points).
    bool lowIsChain = false;
    bool highIsChain = false;
    // Edge of the cut copy between the regions with the shortest way from
    // source to source; -1 where two chains meet in a point.
    int edge = -1;
    // Coordinate the regions were grown for.
    int coordinate = 0;
    // Variables of that coordinate at the sources the two ends of the edge
    // belong to, in the order of the ends (-1 without edge).
    int first = -1;
    int second = -1;
    // Distance of the two sources in the coordinate.
    float span = 0.0f;
};

// The bridges of a part, those of the first coordinate in front of those of
// the second. The number of a bridge is its place in this list.
class RegionBridges {
public:
    // Appends the bridges between the regions of `growth`: every edge of the
    // copy whose ends lie in two regions is a candidate for their pair, in
    // ascending order; the cheaper one wins, the earlier one on ties. An edge
    // of a chain between two sources joins two chains that meet: it always
    // wins and leaves the pair without edge.
    void add(const RegionGrowth& growth, const UvOutput::Coordinates& uv, const mesh::Adjacency& copy);

    [[nodiscard]] int count() const noexcept { return static_cast<int>(list_.size()); }
    [[nodiscard]] const std::vector<Bridge>& list() const noexcept { return list_; }

    // The numbers of all bridges ordered by span, smallest first.
    [[nodiscard]] std::vector<int> bySpan() const;

private:
    std::vector<Bridge> list_;
};

} // namespace areko::remesher::param
