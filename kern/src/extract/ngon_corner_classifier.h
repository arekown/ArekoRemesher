#pragma once

#include "ngon_new_points.h"
#include "ngon_polygon_list.h"

namespace areko::remesher::extract {

// Flag counts of one polygon: corners (1), straights (0), concave (-1),
// spikes (2), unclear (3), and the last position with flag 0 (-1 if none).
struct FlagCounts {
    int corners = 0;
    int straights = 0;
    int concave = 0;
    int spikes = 0;
    int unclear = 0;
    int lastStraight = -1;
};

// Classifies the corners of polygons by their turning angle around the
// polygon normal.
class CornerClassifier {
public:
    // Grows the corner data to the corner count of the list (angle -100, flag 3)
    // and classifies every polygon that ends behind the corner data it had.
    static void classify(const NgonPolygonList& list, const NgonNewPoints& points, NgonCornerData& data);
    [[nodiscard]] static FlagCounts count(const NgonCornerData& data, int first, int last);
    // Second look at unclear corners: with four corners and unclear ones but no
    // concave one, unclear corners below 153 degrees become corners; then, with
    // exactly four corners and no concave one, unclear corners between 135 and
    // 225 degrees become straight.
    static void reclassify(NgonCornerData& data, FlagCounts& counts, int first, int last);

private:
    static void classifyPolygon(const NgonPolygonList& list, const NgonNewPoints& points, int p, NgonCornerData& data);
};

} // namespace areko::remesher::extract
