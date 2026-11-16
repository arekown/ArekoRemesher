#include "curvature_radii.h"

#include "util/select.h"

#include <bit>
#include <cmath>
#include <cstdint>

namespace areko::remesher::field {

namespace {

// Factor of the smallest curvature, about 1 / (2 tan 7.5 degrees), as the float the computation uses.
constexpr float minCurvatureFactor = 3.7978768f;
static_assert(std::bit_cast<std::uint32_t>(minCurvatureFactor) == 0x4073106au);
static_assert(std::bit_cast<std::uint32_t>(0.33f) == 0x3ea8f5c3u && std::bit_cast<std::uint32_t>(1.7f) == 0x3fd9999au);

} // namespace

CurvatureRadii CurvatureRadii::compute(float edgeLength, float area, bool borderData, bool severalParts)
{
    CurvatureRadii r;
    float length = edgeLength;
    if (!borderData) {
        length = util::maxss(std::sqrt(area / 10000.0f), edgeLength);
    }
    r.rMin_ = 0.33f * length;
    r.rMax_ = 1.75f * length;
    if (static_cast<double>(r.rMin_) > static_cast<double>(r.rMax_) * 0.8) {
        r.rMin_ = static_cast<float>(static_cast<double>(r.rMax_) * 0.8);
    }
    const float q = r.rMax_ / r.rMin_;
    if (severalParts) {
        const float c = std::sqrt(area / 24.0f) * 0.5f;
        if (c < r.rMax_) {
            r.rMax_ = c;
        }
        if (c < r.rMin_) {
            r.rMin_ = c;
        }
        if (q == 1.0f) {
            r.levels_ = 1;
            return r;
        }
        if (r.rMax_ == c && c == r.rMin_) {
            r.rMin_ = r.rMax_ / q;
        }
    } else if (q == 1.0f) {
        r.levels_ = 1;
        return r;
    }
    if (!(1.7f <= q)) {
        r.levels_ = 2;
    } else if (q < 3.0f) {
        r.levels_ = 3;
    } else if (!(4.0f <= q)) {
        r.levels_ = 4;
    } else {
        r.levels_ = 7.0f <= q ? 6 : 5;
    }
    return r;
}

float CurvatureRadii::levelRadius(int level) const
{
    if (levels_ == 1) {
        return rMax_;
    }
    const float t = static_cast<float>(level) / static_cast<float>(levels_ - 1);
    const float g = rMax_ / std::pow(10.0f, (1.0f - t) * std::log10(rMax_ / rMin_));
    return (((rMax_ - rMin_) * t + rMin_) + g) * 0.5f;
}

float CurvatureRadii::minCurvature(float radius)
{
    return 1.0f / (radius * minCurvatureFactor);
}

} // namespace areko::remesher::field
