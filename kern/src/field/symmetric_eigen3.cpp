#include "symmetric_eigen3.h"

#include <cmath>
#include <utility>

namespace areko::remesher::field {

namespace {

constexpr double tolerance = 1e-20;
constexpr int maxSweeps = 70;

// One plane rotation applied to the pair (g, h).
struct Rotation {
    double s = 0.0;
    double tau = 0.0;

    void apply(double& g, double& h) const
    {
        const double g0 = g;
        const double h0 = h;
        g = g0 - (g0 * tau + h0) * s;
        h = (g0 - h0 * tau) * s + h0;
    }
};

} // namespace

std::array<SymmetricEigen3::Pair, 3> SymmetricEigen3::solve(const std::array<double, 6>& m)
{
    double a[3][3] = {{m[0], m[1], m[2]}, {m[1], m[3], m[4]}, {m[2], m[4], m[5]}};
    double v[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    double b[3] = {m[0], m[3], m[5]};
    double d[3] = {m[0], m[3], m[5]};
    double z[3] = {0.0, 0.0, 0.0};

    auto offDiagonal = [&a] { return (std::fabs(a[0][1]) + std::fabs(a[0][2])) + std::fabs(a[1][2]); };
    double sm = offDiagonal();
    if (!(tolerance > sm)) {
        for (int sweep = 1;; ++sweep) {
            const double tresh = sweep < 4 ? (sm * 0.2) / 9.0 : 0.0;
            for (int p = 0; p < 2; ++p) {
                for (int q = p + 1; q < 3; ++q) {
                    const double g = std::fabs(a[p][q]) * 100.0;
                    if (sweep > 4 && tolerance > g) {
                        a[p][q] = 0.0;
                        continue;
                    }
                    if (!(std::fabs(a[p][q]) > tresh)) {
                        continue;
                    }
                    const double h = d[q] - d[p];
                    double t;
                    if (tolerance > g) {
                        t = a[p][q] / h;
                    } else {
                        const double theta = (h * 0.5) / a[p][q];
                        t = 1.0 / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                        if (theta < 0.0) {
                            t = -t;
                        }
                    }
                    const double c = 1.0 / std::sqrt(t * t + 1.0);
                    Rotation r;
                    r.s = t * c;
                    r.tau = r.s / (1.0 + c);
                    const double shift = t * a[p][q];
                    a[p][q] = 0.0;
                    z[p] = z[p] - shift;
                    z[q] = z[q] + shift;
                    d[p] = d[p] - shift;
                    d[q] = shift + d[q];
                    // The other two off-diagonal entries as one pair.
                    if (p == 0 && q == 1) {
                        r.apply(a[0][2], a[1][2]);
                    } else if (p == 0) {
                        r.apply(a[0][1], a[1][2]);
                    } else {
                        r.apply(a[0][1], a[0][2]);
                    }
                    for (auto& row : v) {
                        r.apply(row[p], row[q]);
                    }
                }
            }
            for (int i = 0; i < 3; ++i) {
                b[i] = b[i] + z[i];
                d[i] = b[i];
                z[i] = 0.0;
            }
            sm = offDiagonal();
            if (sweep == maxSweeps || tolerance > sm) {
                break;
            }
        }
    }

    std::array<Pair, 3> out;
    for (int k = 0; k < 3; ++k) {
        out[static_cast<std::size_t>(k)] = {{v[0][k], v[1][k], v[2][k]}, d[k]};
    }
    auto order = [&out](std::size_t i, std::size_t j) {
        if (std::fabs(out[j].value) < std::fabs(out[i].value)) {
            std::swap(out[i], out[j]);
        }
    };
    order(0, 1);
    order(1, 2);
    order(0, 1);
    return out;
}

} // namespace areko::remesher::field
