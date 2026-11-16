#pragma once

#include "guide_curve.h"
#include "polyline_locator.h"

#include "mesh/mesh.h"
#include "parts/part_builder.h"
#include "prep/net_lines.h"

#include <span>
#include <vector>

namespace areko::remesher::curves {

// Hands the guide curves and the lines of the preparation to the parts:
//  - every line becomes a curve (index -2, type and flags of the line, its
//    points as point locations; only lines with flags exactly 1 or 2 keep
//    them, the others are located anew),
//  - candidates are the host curves followed by the line curves; inactive
//    ones are left out,
//  - a curve without valid locations is laid onto the prepared mesh point by
//    point, snapping to corner points (PolylineLocator),
//  - the located points go to the part of their face (point locations to the
//    part point at the same position); whenever the part changes, the run
//    collected so far goes to its part when it has at least two points.
class CurveDistribution {
public:
    static void run(const mesh::Mesh& prepared, std::vector<GuideCurve> hostCurves,
                    std::span<const prep::NetLine> lines, parts::PartBuilder::Split& split);

private:
    // Curves of the lines of the preparation.
    [[nodiscard]] static std::vector<GuideCurve> fromLines(std::span<const prep::NetLine> lines);

    // Locations of a curve on the prepared mesh: taken over when valid,
    // otherwise found point by point.
    [[nodiscard]] static LocatedPolyline locate(const GuideCurve& curve, PolylineLocator& locator);
};

} // namespace areko::remesher::curves
