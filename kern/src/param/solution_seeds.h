#pragma once

#include "param/solution.h"
#include "param/system_context.h"
#include "prep/symmetry_cache.h"

#include <array>
#include <span>
#include <vector>

namespace areko::remesher::param {

// The start points the extraction gets besides the singular points, and what
// symmetry does to both lists.
class SolutionSeeds {
public:
    // For a part without any singular point: every border point of the cut
    // copy whose two coordinates are whole numbers adds its point of the part
    // mesh to the seeds, once, in the order of the border points. Runs only
    // while the seed list is empty or the borders are to be frozen.
    static void addBorderPoints(const SystemContext& context, bool freezeBorders, Solution& solution);

    // Appends to `seeds` every border point of the cut copy whose two
    // coordinates `uv` are whole numbers, as its point of the part mesh (of
    // `pointCount` points), in the order of the border points and only when
    // that point is not listed yet. The extraction calls it as well.
    static void appendBorderPoints(const field::CutCopy& copy, int pointCount,
                                   std::span<const std::array<float, 2>> uv, std::vector<int>& seeds);

    // On a symmetric part: singular points and seeds that touch only
    // triangles of the mirrored half are dropped; every other one is kept
    // and followed by its mirror point (same index and mark) unless it lies
    // on the plane or has none.
    static void mirror(const SystemContext& context, const prep::SymmetryCache& pairs, Solution& solution);
};

} // namespace areko::remesher::param
