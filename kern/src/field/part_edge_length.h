#pragma once

#include "parts/part.h"

#include <cstddef>

namespace areko::remesher::field {

// Target edge length h of one part for the field steps.
class PartEdgeLength {
public:
    // Writes part.edgeLength and part.bordersFrozen. With freeze marks h is
    // the mean chord between consecutive kept border points, without the
    // chords over locked border edges, divided by `borderDensityInverse`;
    // it counts when above 0, and only then the borders of the part are
    // frozen. Otherwise: without area (not > 0) h stays; N is the part's quad
    // budget, raised when below 20: to 20 when the quad count is above 400 or
    // there are at least three parts, else to at least 6 (no border edge) or
    // 4; h is the given edge length when above 0, else sqrtf(A / N).
    static void assign(parts::Part& part, int quadCount, std::size_t partCount, float borderDensityInverse = 1.0f);
};

} // namespace areko::remesher::field
