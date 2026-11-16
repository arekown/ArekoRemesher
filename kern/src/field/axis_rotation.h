#pragma once

#include "mesh/vec3.h"

#include <array>
#include <cmath>
#include <cstddef>

namespace areko::remesher::field {

// Rotation about an axis by an angle, kept as the nine float entries of its
// matrix (columns in order). `Rule` names what differs from use to use:
//   Real          number type of the axis length, of the four parts made from
//                 the half angle and of their products, and the type in which
//                 apply() sums a coordinate,
//   shortestAxis  axes shorter than this give no rotation,
//   addsZero      whether apply() adds +0 to every coordinate (-0 becomes +0).
// Whatever `Real` is, sine and cosine of the half angle are taken in double,
// and the products are summed and doubled in double.
template <class Rule>
class AxisRotation {
public:
    using Real = typename Rule::Real;

    // No rotation.
    AxisRotation() = default;

    AxisRotation(mesh::Vec3f axis, float angle)
    {
        const Real axisLength = std::sqrt(static_cast<Real>(mesh::dot(axis, axis)));
        if (static_cast<double>(axisLength) < Rule::shortestAxis) {
            return;
        }
        double sine = 0.0;
        double cosine = 0.0;
        halfAngle(mesh::wide(angle) * 0.5, sine, cosine);
        const Real w = static_cast<Real>(cosine);
        const Real q = static_cast<Real>(sine) / axisLength;
        const Real x = static_cast<Real>(axis.x) * q;
        const Real y = static_cast<Real>(axis.y) * q;
        const Real z = q * static_cast<Real>(axis.z);
        const double ww = static_cast<double>(w * w);
        const double wx = static_cast<double>(w * x);
        const double wy = static_cast<double>(w * y);
        const double wz = static_cast<double>(w * z);
        const double xy = static_cast<double>(x * y);
        const double xz = static_cast<double>(x * z);
        const double yz = static_cast<double>(y * z);
        m_[0] = diagonal(static_cast<double>(x * x) + ww);
        m_[1] = twice(wz + xy);
        m_[2] = twice(xz - wy);
        m_[3] = twice(xy - wz);
        m_[4] = diagonal(static_cast<double>(y * y) + ww);
        m_[5] = twice(wx + yz);
        m_[6] = twice(wy + xz);
        m_[7] = twice(yz - wx);
        m_[8] = diagonal(static_cast<double>(z * z) + ww);
    }

    // The rotated vector: every coordinate summed in `Real`, then rounded.
    [[nodiscard]] mesh::Vec3f apply(mesh::Vec3f v) const
    {
        auto row = [&](std::size_t i) {
            const Real sum = (static_cast<Real>(m_[i]) * static_cast<Real>(v.x) +
                              static_cast<Real>(m_[3 + i]) * static_cast<Real>(v.y)) +
                             static_cast<Real>(m_[6 + i]) * static_cast<Real>(v.z);
            return static_cast<float>(Rule::addsZero ? sum + Real{0} : sum);
        };
        return {row(0), row(1), row(2)};
    }

private:
    // Sine and cosine of the half angle from one call where the C library has it.
    static void halfAngle(double half, double& sine, double& cosine)
    {
#if defined(__GLIBC__)
        ::sincos(half, &sine, &cosine);
#else
        sine = std::sin(half);
        cosine = std::cos(half);
#endif
    }
    static float diagonal(double sum) { return static_cast<float>((sum + sum) - 1.0); }
    static float twice(double sum) { return static_cast<float>(sum + sum); }

    std::array<float, 9> m_{1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f};
};

// A plane (or a direction in it) turned into another plane, in double
// throughout: axes shorter than 1e-7 give no rotation, and a coordinate of -0
// comes out as +0. Used for the frames of the cross field and for the walk
// over the triangle mesh.
struct PlaneTurn {
    using Real = double;
    static constexpr double shortestAxis = 1e-7;
    static constexpr bool addsZero = true;
};

} // namespace areko::remesher::field
