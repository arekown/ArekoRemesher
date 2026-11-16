#pragma once

#include "mesh/vec3.h"

#include <array>

namespace areko::remesher::prep {

// Plane quadric in the order a00 a01 a02 a11 a12 a22 b0 b1 b2 c (double).
struct Quadric {
    std::array<double, 10> q{};

    // Area-weighted plane quadric of a triangle; the weight is rounded to float.
    static Quadric ofTriangle(mesh::Vec3f a, mesh::Vec3f b, mesh::Vec3f c, float weightScale)
    {
        const mesh::Vec3d p0{a.x, a.y, a.z};
        const mesh::TriangleNormal normal = mesh::triangleNormal(a, b, c);
        const mesh::Vec3d& n = normal.unit;
        const double len = normal.length;
        const double d = -((p0.x * n.x + p0.y * n.y) + p0.z * n.z);
        const double w = static_cast<double>(static_cast<float>(len * static_cast<double>(weightScale)));
        Quadric out;
        out.q = {n.x * (w * n.x), n.y * (w * n.x), n.z * (w * n.x), n.y * (n.y * w), n.z * (n.y * w),
                 n.z * (n.z * w), (w * n.x) * d,   (n.y * w) * d,   (n.z * w) * d,   (d * d) * w};
        return out;
    }

    Quadric& operator+=(const Quadric& other)
    {
        for (std::size_t i = 0; i < q.size(); ++i) {
            q[i] += other.q[i];
        }
        return *this;
    }
    friend Quadric operator+(Quadric a, const Quadric& b) { return a += b; }

    // Error at a point, summed in a fixed order.
    [[nodiscard]] double evaluate(mesh::Vec3f p) const
    {
        const double x = p.x;
        const double y = p.y;
        const double z = p.z;
        return (((((((((q[0] * x) * x + ((q[1] + q[1]) * x) * y) + ((q[2] + q[2]) * x) * z) + (q[6] + q[6]) * x) +
                     (q[3] * y) * y) +
                    ((q[4] + q[4]) * y) * z) +
                   (q[7] + q[7]) * y) +
                  (q[5] * z) * z) +
                 (q[8] + q[8]) * z) +
               q[9];
    }
};

} // namespace areko::remesher::prep
