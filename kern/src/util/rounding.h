#pragma once

#include "mesh/vec3.h"

#include <climits>

namespace areko::remesher::util {

// Truncation toward zero like cvttsd2si: NaN and values outside the int range
// give INT_MIN.
[[nodiscard]] inline int truncateToInt(double x) noexcept
{
    if (!(x > static_cast<double>(INT_MIN) - 1.0 && x < static_cast<double>(INT_MAX) + 1.0)) {
        return INT_MIN;
    }
    return static_cast<int>(x);
}

// Rounds half away from zero: (int)(x - 0.5) for x < 0, else (int)(x + 0.5),
// both truncated like truncateToInt. -0.0 and NaN take the + 0.5 branch.
[[nodiscard]] inline int roundHalfAway(double x) noexcept
{
    return truncateToInt(x < 0.0 ? x - 0.5 : x + 0.5);
}

// The same in float: the half is added in float and the sum truncated the
// way the processor does it for float.
[[nodiscard]] inline int roundHalfAway(float x) noexcept
{
    return mesh::truncToInt(x < 0.0f ? x - 0.5f : x + 0.5f);
}

} // namespace areko::remesher::util
