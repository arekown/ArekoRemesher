#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <numbers>

// Small value types and arithmetic helpers. Every helper fixes the order of
// operations and the precision (float or double), because results must be
// reproducible bit for bit.
namespace areko::remesher::mesh {

// pi rounded to float (0x40490fdb).
inline constexpr float piF = std::numbers::pi_v<float>;
static_assert(std::bit_cast<std::uint32_t>(piF) == 0x40490fdbu);
// 2 pi rounded to float (0x40c90fdb).
inline constexpr float fullTurnF = 6.2831855f;
static_assert(std::bit_cast<std::uint32_t>(fullTurnF) == 0x40c90fdbu);
// pi, pi / 2 and 2 pi in double; halving and doubling pi is exact.
inline constexpr double piD = std::numbers::pi;
inline constexpr double halfPiD = 1.5707963267948966;
inline constexpr double fullTurnD = 6.283185307179586;
static_assert(halfPiD + halfPiD == piD && piD + piD == fullTurnD);
static_assert(std::bit_cast<std::uint64_t>(fullTurnD) == 0x401921fb54442d18ull);

// Angle from its cosine in double: pi for <= -1, 0 for >= 1, else acos (ordered
// tests, so a NaN cosine gives NaN). Rounded to float, pi becomes piF.
inline double angleOfCosine(double cosine)
{
    if (cosine <= -1.0) {
        return piD;
    }
    if (cosine >= 1.0) {
        return 0.0;
    }
    return std::acos(cosine);
}

// Angle limits: degrees times (pi / 180) in double; float limits are that value rounded.
constexpr double degD(double degrees) { return degrees * (std::numbers::pi / 180.0); }
constexpr float degF(double degrees) { return static_cast<float>(degD(degrees)); }
// A float value compared or computed in double.
constexpr double wide(float value) { return static_cast<double>(value); }
// An angle in radians as degrees: times (180 / pi) in double, rounded to float.
constexpr float degreesOf(float radians) { return static_cast<float>(wide(radians) * (180.0 / std::numbers::pi)); }

struct Vec3f {
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;

    friend bool operator==(const Vec3f&, const Vec3f&) = default;
};

struct Vec3d {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Color3ub {
    std::uint8_t r = 0, g = 0, b = 0;
    friend bool operator==(const Color3ub&, const Color3ub&) = default;
};

struct Color4ub {
    std::uint8_t r = 0, g = 0, b = 0, a = 0;
    friend bool operator==(const Color4ub&, const Color4ub&) = default;
};

// Coordinate `axis` (0 = x, 1 = y, 2 = z).
inline float& axisValue(Vec3f& v, int axis) { return axis == 0 ? v.x : (axis == 1 ? v.y : v.z); }
inline float axisValue(const Vec3f& v, int axis) { return axis == 0 ? v.x : (axis == 1 ? v.y : v.z); }

inline Vec3f operator-(Vec3f a, Vec3f b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3f operator+(Vec3f a, Vec3f b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3f operator*(Vec3f a, float s) { return {a.x * s, a.y * s, a.z * s}; }

// (ax*bx + ay*by) + az*bz in float.
inline float dot(Vec3f a, Vec3f b) { return (a.x * b.x + a.y * b.y) + a.z * b.z; }

// sqrt((dx*dx + dy*dy) + dz*dz) in float.
inline float length(Vec3f v) { return std::sqrt((v.x * v.x + v.y * v.y) + v.z * v.z); }

inline float distance(Vec3f a, Vec3f b) { return length(a - b); }

// Length of the edge a-b as the mesh measures keep it: differences, squares,
// (dx*dx + dy*dy) + dz*dz and the root in double, rounded to float. The same
// for b-a; a NaN coordinate gives NaN.
inline float edgeLength(Vec3f a, Vec3f b)
{
    const double x = static_cast<double>(a.x) - static_cast<double>(b.x);
    const double y = static_cast<double>(a.y) - static_cast<double>(b.y);
    const double z = static_cast<double>(a.z) - static_cast<double>(b.z);
    return static_cast<float>(std::sqrt((x * x + y * y) + z * z));
}

// t without its part along n: t - n * (t . n), all in float.
inline Vec3f projectOntoPlane(Vec3f t, Vec3f n) { return t - n * dot(t, n); }

// Order by x, then y, then z; a NaN coordinate never counts as smaller or larger.
inline bool lessXyz(Vec3f p, Vec3f q)
{
    if (q.x > p.x) {
        return true;
    }
    if (p.x > q.x) {
        return false;
    }
    if (q.y > p.y) {
        return true;
    }
    if (p.y > q.y) {
        return false;
    }
    return q.z > p.z;
}

// sqrt in double of the float square sum (x*x + y*y) + z*z.
inline double lengthAsDouble(Vec3f v) { return std::sqrt(static_cast<double>((v.x * v.x + v.y * v.y) + v.z * v.z)); }

// Nearest point of the segment a-b to p, all in float:
// t = clamp((p - a).(b - a) / |b - a|^2, 0, 1) (0 for a segment of length 0),
// point = a + (b - a) * t.
struct SegmentPoint {
    Vec3f point;
    float t = 0.0f;
};
inline SegmentPoint nearestOnSegment(Vec3f p, Vec3f a, Vec3f b)
{
    const Vec3f ab = b - a;
    const float len2 = dot(ab, ab);
    const float t = len2 > 0.0f ? std::clamp(dot(p - a, ab) / len2, 0.0f, 1.0f) : 0.0f;
    return {a + ab * t, t};
}

// Closest point of the segment a-b to p and its squared distance, all in
// float: p not beyond a along b - a gives a; p before b (also when that test
// is NaN) gives the foot on the line, with b - a normalized unless its length
// lies within +-1e-14; otherwise b.
struct SegmentClosest {
    Vec3f point;
    float distanceSq = 0.0f;
};
inline SegmentClosest closestOnSegment(Vec3f p, Vec3f a, Vec3f b)
{
    const Vec3f ab = b - a;
    const Vec3f ap = p - a;
    Vec3f q = a;
    if (!(dot(ap, ab) <= 0.0f)) {
        if (!(dot(p - b, ab) >= 0.0f)) {
            const float length = std::sqrt(dot(ab, ab));
            Vec3f unit{};
            if (!(length > -1e-14f && length < 1e-14f)) {
                unit = {ab.x / length, ab.y / length, ab.z / length};
            }
            const float along = dot(ap, unit);
            const Vec3f off = ap - unit * along;
            return {p - off, dot(off, off)};
        }
        q = b;
    }
    const Vec3f d = p - q;
    return {q, dot(d, d)};
}

// A float vector in double.
inline Vec3d wide(Vec3f v) { return {wide(v.x), wide(v.y), wide(v.z)}; }

inline Vec3d operator-(const Vec3d& a, const Vec3d& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }

// (ax*bx + ay*by) + az*bz in double.
inline double dot(const Vec3d& a, const Vec3d& b) { return (a.x * b.x + a.y * b.y) + a.z * b.z; }

// sqrt((x*x + y*y) + z*z) in double.
inline double length(const Vec3d& v) { return std::sqrt(dot(v, v)); }

// a x b = (ay*bz - az*by, az*bx - ax*bz, ax*by - ay*bx).
inline Vec3f cross(Vec3f a, Vec3f b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }

inline Vec3d cross(const Vec3d& a, const Vec3d& b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

// Coordinate axis along the smallest component of n (x and y win ties in this order).
inline Vec3f smallestAxis(Vec3f n)
{
    const float ax = std::fabs(n.x);
    const float ay = std::fabs(n.y);
    const float az = std::fabs(n.z);
    if (ax <= ay && ax <= az) {
        return {1.0f, 0.0f, 0.0f};
    }
    if (ax > ay && ay < az) {
        return {0.0f, 1.0f, 0.0f};
    }
    return {0.0f, 0.0f, 1.0f};
}

// Float normalization: squared length in float, root and division in double.
// Vectors shorter than 1e-28 stay unchanged. Returns the length as float.
inline float normalizeWithLength(Vec3f& v)
{
    const float q = (v.x * v.x + v.y * v.y) + v.z * v.z;
    const double l = std::sqrt(static_cast<double>(q));
    if (!(-1e-28 < l && l < 1e-28)) {
        v.x = static_cast<float>(static_cast<double>(v.x) / l);
        v.y = static_cast<float>(static_cast<double>(v.y) / l);
        v.z = static_cast<float>(static_cast<double>(v.z) / l);
    }
    return static_cast<float>(l);
}

inline void normalizeInPlace(Vec3f& v) { (void)normalizeWithLength(v); }

// Double normalization; vectors shorter than 1e-100 stay unchanged.
// Returns the length before normalizing.
inline double normalizeInPlace(Vec3d& v)
{
    const double len = std::sqrt((v.x * v.x + v.y * v.y) + v.z * v.z);
    if (!(-1e-100 < len && len < 1e-100)) {
        v.x = v.x / len;
        v.y = v.y / len;
        v.z = v.z / len;
    }
    return len;
}

// Normal of the triangle (a, b, c) in double: edges b - a and c - a, their
// cross product, normalized. `length` is the length before normalizing (the
// doubled area).
struct TriangleNormal {
    Vec3d unit;
    double length = 0.0;
};

inline TriangleNormal triangleNormal(Vec3f a, Vec3f b, Vec3f c)
{
    const Vec3d p0{a.x, a.y, a.z};
    const Vec3d e1{static_cast<double>(b.x) - p0.x, static_cast<double>(b.y) - p0.y, static_cast<double>(b.z) - p0.z};
    const Vec3d e2{static_cast<double>(c.x) - p0.x, static_cast<double>(c.y) - p0.y, static_cast<double>(c.z) - p0.z};
    TriangleNormal out;
    out.unit = cross(e1, e2);
    out.length = normalizeInPlace(out.unit);
    return out;
}

// Float to int32 the way the processor truncates (cvttss2si): towards zero,
// NaN and out-of-range values give INT32_MIN.
inline std::int32_t truncToInt(float value)
{
    if (!(value > -2147483904.0f && value < 2147483648.0f)) {
        return INT32_MIN;
    }
    return static_cast<std::int32_t>(value);
}

// Whether a float is a whole number: truncated and turned back it is the
// same value (false for NaN, infinities and values outside the int range;
// rounding half away from zero instead of truncating decides the same).
inline bool isWholeNumber(float value) { return value == static_cast<float>(truncToInt(value)); }

inline std::int64_t truncToInt64(float value)
{
    if (!(value > -9223373136366403584.0f && value < 9223372036854775808.0f)) {
        return INT64_MIN;
    }
    return static_cast<std::int64_t>(value);
}

inline std::int32_t truncToInt(double value)
{
    if (!(value > -2147483649.0 && value < 2147483648.0)) {
        return INT32_MIN;
    }
    return static_cast<std::int32_t>(value);
}

} // namespace areko::remesher::mesh
