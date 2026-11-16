#pragma once

namespace areko::remesher::util {

// The smaller resp. larger of two values with a fixed choice for ties and
// values that are not a number: a only when it is strictly smaller resp.
// larger than b, in every other case b.
[[nodiscard]] constexpr float minss(float a, float b) noexcept
{
    return a < b ? a : b;
}
[[nodiscard]] constexpr float maxss(float a, float b) noexcept
{
    return a > b ? a : b;
}
[[nodiscard]] constexpr double maxsd(double a, double b) noexcept
{
    return a > b ? a : b;
}

} // namespace areko::remesher::util
