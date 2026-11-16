#include "ngon_cut_check.h"

#include <bit>
#include <cmath>
#include <cstdint>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

constexpr float parallelLimit = std::bit_cast<float>(0x3727c5acu); // 1e-5
constexpr float segmentLow = std::bit_cast<float>(0x38d1b717u);    // 1e-4
constexpr float segmentHigh = std::bit_cast<float>(0x3f7ff972u);   // 0.9999
constexpr double nearFactor = 0.4;
constexpr double behindFactor = 0.05;
constexpr double beyondFactor = 1.05;
constexpr float quarterTurn = std::bit_cast<float>(0x3f490fdbu);   // pi / 4
constexpr float secondLimit = std::bit_cast<float>(0x283424dcu);   // 1e-14

// Angle between the line along u and d: 0 unless |u . d| < 1, else acos (also NaN).
float lineTurn(Vec3f u, Vec3f d)
{
    double x = static_cast<double>(dot(u, d));
    if (0.0 > x) {
        x = -x;
    }
    return x >= 1.0 ? 0.0f : static_cast<float>(std::acos(x));
}

float endFactor(Vec3f before, Vec3f here, Vec3f after, Vec3f d, float cost)
{
    Vec3f a1 = here - before;
    Vec3f a2 = here - after;
    normalizeInPlace(a1);
    normalizeInPlace(a2);
    const float spread = std::fabs(lineTurn(a1, d) - lineTurn(a2, d));
    return (((spread * 0.5f) * 3.0f) / quarterTurn + 1.0f) * cost;
}

} // namespace

Approach CutCheck::closest(Vec3f o, Vec3f d, Vec3f a, Vec3f b)
{
    const Vec3f edge = b - a;
    const float span = length(edge);
    const Vec3f e{edge.x / span, edge.y / span, edge.z / span};
    const Vec3f c = cross(d, e);
    const float cl = length(c);
    const Vec3f w = a - o;
    Approach out;
    float s = 0.0f;
    bool parallel = -parallelLimit < cl && cl < parallelLimit;
    if (!parallel) {
        const float bd = dot(e, d);
        const float dw = dot(d, w);
        const float ew = dot(e, w);
        const float den = -1.0f - bd * -bd;
        if (den != 0.0f) {
            s = (ew - bd * dw) / den;
            out.along = (-dw - -bd * ew) / den;
            out.distance = std::fabs(dot(w, c) / cl);
        } else {
            parallel = true;
        }
    }
    if (parallel) {
        s = 0.0f;
        out.along = dot(d, w);
        out.distance = length(cross(d, w));
    }
    float u = s / span;
    Vec3f end;
    if (u < 0.0f) {
        u = 0.0f;
        end = a;
    } else if (!(u > 1.0f)) {
        // Also a parameter that is not a number (segment of length zero).
        out.segment = u;
        return out;
    } else {
        u = 1.0f;
        end = b;
    }
    out.segment = u;
    const Vec3f q = end - o;
    out.along = dot(d, q);
    out.distance = length(q - d * out.along);
    return out;
}

CutCheck::Chord CutCheck::chord(const NgonPolygonList& list, const NgonNewPoints& points, int p, int i, int j)
{
    const Vec3f from = points.position(list.vertex(i));
    Vec3f d = points.position(list.vertex(j)) - from;
    Chord out;
    out.length = normalizeWithLength(d);
    // The direction is normalized a second time in float: by the root of its
    // squared length, (0, 0, 0) only while that root lies strictly between
    // -1e-14 and 1e-14 (a root that is not a number divides).
    const float again = length(d);
    if (-secondLimit < again && again < secondLimit) {
        d = {};
    } else {
        d = {d.x / again, d.y / again, d.z / again};
    }
    const double limit = static_cast<double>(out.length);
    for (int r = list.begin(p); r < list.end(p); ++r) {
        const int r1 = list.next(p, r);
        if (r == i || r1 == i || r == j || r1 == j) {
            continue;
        }
        const Approach near =
            closest(from, d, points.position(list.vertex(r)), points.position(list.vertex(r1)));
        if (near.segment > segmentLow && near.segment < segmentHigh &&
            static_cast<double>(near.distance) < nearFactor * limit && -behindFactor * limit < static_cast<double>(near.along) &&
            static_cast<double>(near.along) < beyondFactor * limit) {
            out.blocked = true;
            return out;
        }
    }
    return out;
}

float CutCheck::angleFactor(const NgonPolygonList& list, const NgonNewPoints& points, int p, int i, int j, float length,
                            float cost)
{
    const Vec3f pi = points.position(list.vertex(i));
    const Vec3f pj = points.position(list.vertex(j));
    const Vec3f diff = pi - pj;
    const Vec3f d{diff.x / length, diff.y / length, diff.z / length};
    cost = endFactor(points.position(list.vertex(list.prev(p, i))), pi, points.position(list.vertex(list.next(p, i))),
                     d, cost);
    return endFactor(points.position(list.vertex(list.prev(p, j))), pj, points.position(list.vertex(list.next(p, j))), d,
                     cost);
}

} // namespace areko::remesher::extract
