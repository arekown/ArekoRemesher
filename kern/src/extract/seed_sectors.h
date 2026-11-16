#pragma once

#include "trace_frame.h"

#include <array>

namespace areko::remesher::extract {

// Where a face continues at a seed node on a point of the mesh: the
// triangles around the point are visited from the one the trace arrived in,
// and the first sector the grid line can leave in a growing direction is
// taken.
class SeedSectors {
public:
    // `triangle` (in: arrival, out: chosen; -1 when the fan ends at a missing
    // edge), `direction` (in: arrival, out: 0 or 1, or the arrival turned by
    // a quarter when the fan closes without choice) and `uv` (out: grid
    // coordinates of the point in the chosen triangle, not written when the
    // choice comes from the round). False when the point is not found or has
    // no whole coordinates, or when the round finds nothing.
    static bool choose(const TraceFrame& frame, int point, int& triangle, int& direction, std::array<float, 2>& uv);
};

} // namespace areko::remesher::extract
