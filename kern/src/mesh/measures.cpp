#include "measures.h"

#include "mesh.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>

namespace areko::remesher::mesh {

namespace {

// Cosine of the angle between two vectors in double: a . b / (|a| |b|), none
// when |a| |b| is 0.
std::optional<double> cosineBetween(Vec3f a, Vec3f b)
{
    const double denominator = lengthAsDouble(a) * lengthAsDouble(b);
    if (denominator == 0.0) {
        return std::nullopt;
    }
    return static_cast<double>(dot(a, b)) / denominator;
}

// Angle between two vectors in double: zero when |a| |b| is 0, else acos of
// the quotient, pi for <= -1 and 0 for >= 1.
double angleBetween(Vec3f a, Vec3f b)
{
    const std::optional<double> cosine = cosineBetween(a, b);
    if (!cosine) {
        return 0.0;
    }
    return angleOfCosine(*cosine);
}

} // namespace

std::vector<float> cornerAngles(const Mesh& mesh)
{
    const Topology& topology = mesh.topology();
    if (topology.mode() != FaceMode::Triangles) {
        return {};
    }
    const auto& points = mesh.points();
    const auto& corners = topology.cornerVertex();
    std::vector<float> angles(corners.size(), 0.0f);
    for (std::size_t c = 0; c + 2 < corners.size(); c += 3) {
        const Vec3f p0 = points[static_cast<std::size_t>(corners[c])];
        const Vec3f p1 = points[static_cast<std::size_t>(corners[c + 1])];
        const Vec3f p2 = points[static_cast<std::size_t>(corners[c + 2])];
        const Vec3f u = p2 - p0;
        const Vec3f w = p1 - p0;
        const Vec3f z = p1 - p2;
        const double a0 = angleBetween(u, w);
        const double a1 = angleBetween(w, z);
        angles[c] = static_cast<float>(a0);
        angles[c + 1] = static_cast<float>(a1);
        angles[c + 2] = static_cast<float>((std::numbers::pi - a0) - a1);
    }
    return angles;
}

float vectorAngle(Vec3f u, Vec3f v)
{
    return static_cast<float>(angleBetween(u, v));
}

float lineAngle(Vec3f u, Vec3f v)
{
    const std::optional<double> cosine = cosineBetween(u, v);
    if (!cosine) {
        return 0.0f;
    }
    double quotient = *cosine;
    if (quotient < 0.0) {
        quotient = -quotient;
    }
    return quotient >= 1.0 ? 0.0f : static_cast<float>(std::acos(quotient));
}

double signedVectorAngleAsDouble(Vec3f u, Vec3f v, Vec3f crossed, Vec3f normal)
{
    const double a = angleBetween(u, v);
    if (dot(crossed, normal) > 0.0f) {
        return a;
    }
    // pi stays pi: 2 pi - pi is exact.
    return a != 0.0 ? fullTurnD - a : 0.0;
}

float signedVectorAngle(Vec3f u, Vec3f v, Vec3f crossed, Vec3f normal)
{
    return static_cast<float>(signedVectorAngleAsDouble(u, v, crossed, normal));
}

float cornerAngle(const Mesh& mesh, int corner)
{
    const auto& points = mesh.points();
    const Topology& topology = mesh.topology();
    const int next = corner + (corner % 3 == 2 ? -2 : 1);
    const int prev = corner + (corner % 3 == 0 ? 2 : -1);
    const Vec3f p = points[static_cast<std::size_t>(topology.vertex(corner))];
    return vectorAngle(points[static_cast<std::size_t>(topology.vertex(prev))] - p,
                       points[static_cast<std::size_t>(topology.vertex(next))] - p);
}

std::vector<float> angleDefects(const Adjacency& adjacency, const std::vector<float>& cornerAngles)
{
    const int n = adjacency.pointCount();
    std::vector<float> defects(static_cast<std::size_t>(n), fullTurnF);
    for (int v = 0; v < n; ++v) {
        const auto corners = adjacency.pointCorners(v);
        if (corners.empty()) {
            continue;
        }
        // The float angles are summed in double and subtracted from 2 pi in
        // double; only the difference is rounded to float.
        double sum = 0.0;
        for (const int c : corners) {
            sum += static_cast<double>(cornerAngles[static_cast<std::size_t>(c)]);
        }
        defects[static_cast<std::size_t>(v)] = static_cast<float>(fullTurnD - sum);
    }
    for (const int v : adjacency.boundaryPoints()) {
        defects[static_cast<std::size_t>(v)] = 0.0f;
    }
    return defects;
}

std::vector<float> directAngleDefects(const Mesh& mesh)
{
    const Adjacency& adjacency = mesh.adjacency();
    const auto& points = mesh.points();
    const auto& corners = mesh.topology().cornerVertex();
    const int n = adjacency.pointCount();
    std::vector<float> defects(static_cast<std::size_t>(n), fullTurnF);
    for (int v = 0; v < n; ++v) {
        const auto own = adjacency.pointCorners(v);
        if (own.empty()) {
            continue;
        }
        double sum = 0.0;
        for (const int c : own) {
            const int first = c - c % 3;
            const Vec3f p = points[static_cast<std::size_t>(corners[static_cast<std::size_t>(c)])];
            const Vec3f next = points[static_cast<std::size_t>(corners[static_cast<std::size_t>(first + (c - first + 1) % 3)])];
            const Vec3f previous = points[static_cast<std::size_t>(corners[static_cast<std::size_t>(first + (c - first + 2) % 3)])];
            sum += angleBetween(next - p, previous - p);
        }
        defects[static_cast<std::size_t>(v)] = static_cast<float>(fullTurnD - sum);
    }
    for (const int v : adjacency.boundaryPoints()) {
        defects[static_cast<std::size_t>(v)] = 0.0f;
    }
    return defects;
}

std::vector<float> edgeLengths(const Mesh& mesh)
{
    const Adjacency& adjacency = mesh.adjacency();
    const auto& points = mesh.points();
    std::vector<float> lengths(static_cast<std::size_t>(adjacency.edgeCount()));
    for (int e = 0; e < adjacency.edgeCount(); ++e) {
        const auto ends = adjacency.edgeVertices(e);
        lengths[static_cast<std::size_t>(e)] =
            edgeLength(points[static_cast<std::size_t>(ends[0])], points[static_cast<std::size_t>(ends[1])]);
    }
    return lengths;
}

void addEdgeEndPoints(const Adjacency& adjacency, const BitArray& edges, int pointCount, BitArray& points)
{
    points.resize(static_cast<std::size_t>(pointCount));
    const std::size_t n = std::min(edges.size(), static_cast<std::size_t>(adjacency.edgeCount()));
    for (std::size_t e = 0; e < n; ++e) {
        if (edges.test(e)) {
            const auto ends = adjacency.edgeVertices(static_cast<int>(e));
            points.set(static_cast<std::size_t>(ends[0]));
            points.set(static_cast<std::size_t>(ends[1]));
        }
    }
}

} // namespace areko::remesher::mesh
