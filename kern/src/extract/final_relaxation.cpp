#include "final_relaxation.h"

#include "surface_walker.h"
#include "trace_frame.h"
#include "trace_node.h"

#include "mesh/normals.h"
#include "mesh/vertex_fan.h"
#include "parts/part.h"
#include "prep/polygon_triangulator.h"

#include <cmath>
#include <limits>
#include <optional>
#include <span>

namespace areko::remesher::extract {

using mesh::BitArray;
using mesh::Vec3d;
using mesh::Vec3f;

namespace {

constexpr int rounds = 3;
// Share of the way to the target per step: on a line or a regular point, and for the mean of the neighbours.
constexpr float pullFactor = 0.33f;
constexpr float meanFactor = 0.4f;
// Two dihedral angles above this sum (30 degrees) keep the cross pair of a line point out.
constexpr double creaseLimit = 0.5235987755982988;
// Anisotropy below the first limit counts as none, above the second as full.
constexpr double anisotropyLow = 0.2;
constexpr double anisotropyHigh = 0.85;
// Two fixed edges at a point that bend more than this (cosine) fix the point.
constexpr double bendLimit = 0.7;
constexpr int triangleSource = 4;

// One end of an edge: its length range, its weight, and its direction with that direction's weight.
struct EdgeEnd {
    float lengthMin = 0.0f;
    float lengthMax = 0.0f;
    float weight = 0.0f;
    Vec3f direction;
    float directionWeight = 0.0f;
};

// Value of entry i of a list that may be shorter than the points (-1 beyond).
int bounded(const std::vector<int>& list, int i)
{
    return i >= 0 && static_cast<std::size_t>(i) < list.size() ? list[static_cast<std::size_t>(i)] : -1;
}

// The edge targets of the quad mesh, read from the fields of the part.
class EdgeTargetBuilder {
public:
    EdgeTargetBuilder(parts::Part& part, const RawMesh& raw)
        : part_(part), raw_(raw), triangles_(part.mesh.adjacency()), triPoints_(part.mesh.points()),
          points_(raw.mesh.points()), frame_(part, *part.trace)
    {
    }

    [[nodiscard]] FinalRelaxation::EdgeTargets build()
    {
        const mesh::Adjacency& adjacency = raw_.mesh.adjacency();
        const int edges = adjacency.edgeCount();
        FinalRelaxation::EdgeTargets out;
        out.direction.assign(static_cast<std::size_t>(edges), Vec3f{});
        out.lengthMin.assign(static_cast<std::size_t>(edges), 0.0f);
        out.lengthMax.assign(static_cast<std::size_t>(edges), 0.0f);
        for (int e = 0; e < edges; ++e) {
            const auto ends = adjacency.edgeVertices(e);
            const EdgeEnd p = end(ends[0], ends[1]);
            const EdgeEnd q = end(ends[1], ends[0]);
            const auto i = static_cast<std::size_t>(e);
            const float f = p.weight / (q.weight + p.weight);
            out.lengthMin[i] = p.lengthMin * f + q.lengthMin * (1.0f - f);
            out.lengthMax[i] = f * p.lengthMax + (1.0f - f) * q.lengthMax;
            if (q.directionWeight + p.directionWeight != 0.0f) {
                const float g = p.directionWeight / (q.directionWeight + p.directionWeight);
                const float h = 1.0f - g;
                Vec3f r{p.direction.x * g - q.direction.x * h, p.direction.y * g - q.direction.y * h,
                        p.direction.z * g - q.direction.z * h};
                mesh::normalizeInPlace(r);
                out.direction[i] = r;
            }
            if (raw_.curveEdges.testBounded(i)) {
                Vec3f r = points_[static_cast<std::size_t>(ends[1])] - points_[static_cast<std::size_t>(ends[0])];
                mesh::normalizeInPlace(r);
                out.direction[i] = r;
            }
        }
        return out;
    }

private:
    [[nodiscard]] PointSource sourceOf(int p) const
    {
        const std::vector<PointSource>& s = raw_.points.source;
        return static_cast<std::size_t>(p) < s.size() ? s[static_cast<std::size_t>(p)] : PointSource{};
    }
    [[nodiscard]] std::uint8_t flags(int t) const
    {
        return part_.curvature ? part_.curvature->flags[static_cast<std::size_t>(t)] : std::uint8_t{0};
    }

    // Lengths from the sizes of the triangle of the source, else from the target edge length.
    [[nodiscard]] EdgeEnd fallback(PointSource s) const
    {
        // A point or edge without corners has no triangle.
        auto firstTriangle = [](std::span<const int> corners) { return corners.empty() ? -1 : corners[0] / 3; };
        int t = -1;
        if (s.type == node_kind::point) {
            t = firstTriangle(triangles_.pointCorners(s.id));
        } else if (s.type == node_kind::edge) {
            t = firstTriangle(triangles_.edgeCorners(s.id));
        } else if (s.type == triangleSource) {
            t = s.id;
        }
        EdgeEnd out;
        if (t >= 0) {
            const float su = part_.sizeField->sizeU[static_cast<std::size_t>(t)];
            const float sv = part_.sizeField->sizeV[static_cast<std::size_t>(t)];
            const float mean = (su + sv) * 0.5f;
            out.lengthMin = mean * 0.7f;
            out.lengthMax = mean * 1.5f;
            if (su < out.lengthMin) {
                out.lengthMin = su * 0.95f;
            }
            if (sv < out.lengthMin) {
                out.lengthMin = sv * 0.95f;
            }
            if (out.lengthMax < su) {
                out.lengthMax = su * 1.1f;
            }
            if (out.lengthMax < sv) {
                out.lengthMax = sv * 1.1f;
            }
        } else {
            out.lengthMin = part_.edgeLength * 0.6f;
            out.lengthMax = part_.edgeLength * 1.8f;
        }
        out.weight = 0.2f;
        return out;
    }

    // The end at p of the edge p-q.
    [[nodiscard]] EdgeEnd end(int p, int q)
    {
        const PointSource source = sourceOf(p);
        int t = bounded(raw_.points.triangle, p);
        if (t < 0) {
            return fallback(source);
        }
        int slot = -1;
        for (int k = 0; k < 4 && slot < 0; ++k) {
            if (bounded(raw_.points.neighbours[static_cast<std::size_t>(k)], p) == q) {
                slot = k;
            }
        }
        if (slot < 0) {
            return fallback(source);
        }
        bool flipped = frame_.flipped(t);
        TraceNode node;
        node.kind = source.type;
        node.reference = source.id;
        const Vec3f towards = points_[static_cast<std::size_t>(q)] - points_[static_cast<std::size_t>(p)];
        int better = -1; // a triangle of the source that lies better along the edge
        if (source.type == node_kind::edge) {
            const int e = source.id;
            int c = -1;
            for (const int corner : triangles_.edgeCorners(e)) {
                if (corner / 3 == t) {
                    c = corner;
                    break;
                }
            }
            const int other = c >= 0 ? triangles_.otherCorner(e, c) : -1;
            if (other >= 0) {
                const auto ends = triangles_.edgeVertices(e);
                const Vec3f middle = (triPoints_[static_cast<std::size_t>(ends[0])] +
                                      triPoints_[static_cast<std::size_t>(ends[1])]) *
                                     0.5f;
                const Vec3f own = triPoints_[static_cast<std::size_t>(triangles_.vertex(triangles_.prev(c)))];
                const Vec3f far = triPoints_[static_cast<std::size_t>(triangles_.vertex(triangles_.prev(other)))];
                if (mesh::dot(own - middle, towards) < mesh::dot(far - middle, towards)) {
                    better = other / 3;
                }
            }
        } else if (source.type == node_kind::point) {
            const Vec3f at = triPoints_[static_cast<std::size_t>(source.id)];
            float best = -std::numeric_limits<float>::max();
            for (const int c : triangles_.pointCorners(source.id)) {
                const Vec3f a = triPoints_[static_cast<std::size_t>(triangles_.vertex(triangles_.next(c)))];
                const Vec3f b = triPoints_[static_cast<std::size_t>(triangles_.vertex(triangles_.prev(c)))];
                const float value = mesh::dot((a + b) * 0.5f - at, towards);
                if (best < value) {
                    best = value;
                    better = c / 3;
                }
            }
            if (better == t) {
                better = -1;
            }
        }
        if (better >= 0) {
            const std::optional<int> turned = frame_.directionSlot(better, node, slot, t);
            if (!turned) {
                return fallback(source);
            }
            t = better;
            flipped = flipped || frame_.flipped(t);
            slot = *turned;
            if (slot < 0) {
                return fallback(source);
            }
        }
        const auto i = static_cast<std::size_t>(t);
        const bool second = (slot & 1) != 0;
        const float size = second ? part_.sizeField->sizeV[i] : part_.sizeField->sizeU[i];
        const std::array<float, 2>& c = second ? part_.cutGraph->combedV[i] : part_.cutGraph->combedU[i];
        const Vec3f e1 = part_.crossField->frameE1[i];
        const Vec3f e2 = part_.crossField->frameE2[i];
        EdgeEnd out;
        if (slot < 2) {
            out.direction = {c[1] * e2.x + c[0] * e1.x, c[1] * e2.y + c[0] * e1.y, c[1] * e2.z + c[0] * e1.z};
        } else {
            const float n = -c[0];
            out.direction = {n * e1.x - c[1] * e2.x, n * e1.y - c[1] * e2.y, n * e1.z - c[1] * e2.z};
        }
        out.directionWeight = 1.0f;
        if (flipped && mesh::dot(towards, out.direction) < 0.0f) {
            out.direction = {-out.direction.x, -out.direction.y, -out.direction.z};
        }
        out.lengthMin = size * 0.5f;
        out.lengthMax = (flags(t) & 1u) != 0 ? size * 1.2f : size * 2.5f;
        out.weight = 1.0f;
        return out;
    }

    parts::Part& part_;
    const RawMesh& raw_;
    const mesh::Adjacency& triangles_;
    const std::vector<Vec3f>& triPoints_;
    const std::vector<Vec3f>& points_;
    TraceFrame frame_;
};

// (|l2| - |l1|) / |l2| of triangle t, 0 for |l2| = 0.
float anisotropyOf(const field::Curvature& curvature, int t)
{
    const float k2 = std::fabs(curvature.lambda2[static_cast<std::size_t>(t)]);
    if (k2 == 0.0f) {
        return 0.0f;
    }
    return (k2 - std::fabs(curvature.lambda1[static_cast<std::size_t>(t)])) / k2;
}

// Anisotropy per output point from the curvature around its source.
std::vector<float> pointAnisotropy(const parts::Part& part, const RawMesh& raw)
{
    const auto count = static_cast<std::size_t>(raw.mesh.pointCount());
    std::vector<float> out(count, 0.0f);
    if (!part.curvature) {
        return out;
    }
    const field::Curvature& curvature = *part.curvature;
    const mesh::Adjacency& triangles = part.mesh.adjacency();
    auto marked = [&](int t) { return (curvature.flags[static_cast<std::size_t>(t)] & 1u) != 0; };
    for (std::size_t p = 0; p < count && p < raw.points.source.size(); ++p) {
        const PointSource s = raw.points.source[p];
        float w = 0.0f;
        if (s.type == triangleSource) {
            if (!marked(s.id)) {
                continue;
            }
            w = anisotropyOf(curvature, s.id);
        } else if (s.type == node_kind::point || s.type == node_kind::edge) {
            const auto corners = s.type == node_kind::point ? triangles.pointCorners(s.id) : triangles.edgeCorners(s.id);
            float sum = 0.0f;
            for (const int c : corners) {
                if (marked(c / 3)) {
                    sum = anisotropyOf(curvature, c / 3) + sum;
                }
            }
            w = sum / static_cast<float>(corners.size());
        } else {
            continue;
        }
        if (static_cast<double>(w) < anisotropyLow) {
            w = 0.0f;
        } else if (static_cast<double>(w) > anisotropyHigh) {
            w = 1.0f;
        }
        out[p] = w;
    }
    return out;
}

// Edges of the triangle mesh the walk keeps to: curve edges, border edges and,
// with the symmetry on a symmetric part, the edges on the mirror plane.
BitArray fixedEdges(parts::Part& part, bool symmetry)
{
    BitArray fixed = part.trace->curves.edges;
    const mesh::Adjacency& adjacency = part.mesh.adjacency();
    const auto edgeCount = static_cast<std::size_t>(adjacency.edgeCount());
    const std::vector<int>& border = adjacency.boundaryEdges();
    if (!border.empty()) {
        fixed.resize(edgeCount);
        for (const int e : border) {
            fixed.set(static_cast<std::size_t>(e));
        }
    }
    if (symmetry) {
        if (prep::SymmetryCache* mirror = part.mirror.symmetric(part.mesh, true)) {
            fixed.resize(edgeCount);
            for (const int e : mirror->planeEdges(adjacency)) {
                fixed.set(static_cast<std::size_t>(e));
            }
        }
    }
    return fixed;
}

// 1 per output point, 0 for a kept border point and for a point on a corner of the fixed edges.
std::vector<float> pointMobility(const parts::Part& part, const RawMesh& raw, const BitArray& fixed)
{
    const auto count = static_cast<std::size_t>(raw.mesh.pointCount());
    std::vector<float> out(count, 1.0f);
    const std::vector<PointSource>& sources = raw.points.source;
    const std::size_t known = std::min(count, sources.size());
    if (const std::vector<mesh::FreezeMark>* marks = part.frozenBorderMarks()) {
        for (std::size_t p = 0; p < known; ++p) {
            if (sources[p].type == node_kind::point && (*marks)[static_cast<std::size_t>(sources[p].id)].value == 1) {
                out[p] = 0.0f;
            }
        }
    }
    if (fixed.empty()) {
        return out;
    }
    const mesh::Adjacency& adjacency = part.mesh.adjacency();
    const std::vector<Vec3f>& at = part.mesh.points();
    std::vector<int> edges;
    for (std::size_t p = 0; p < known; ++p) {
        if (sources[p].type != node_kind::point) {
            continue;
        }
        const int v = sources[p].id;
        edges.clear();
        for (const int e : adjacency.pointEdges(v)) {
            if (fixed.testBounded(static_cast<std::size_t>(e))) {
                edges.push_back(e);
            }
        }
        if (edges.size() > 2) {
            out[p] = 0.0f;
        } else if (edges.size() == 2) {
            const int a = adjacency.otherEnd(edges[0], v);
            const int b = adjacency.otherEnd(edges[1], v);
            if (a >= 0 && b >= 0) {
                Vec3f in = at[static_cast<std::size_t>(a)] - at[static_cast<std::size_t>(v)];
                Vec3f outward = at[static_cast<std::size_t>(v)] - at[static_cast<std::size_t>(b)];
                mesh::normalizeInPlace(in);
                mesh::normalizeInPlace(outward);
                if (static_cast<double>(mesh::dot(in, outward)) < bendLimit) {
                    out[p] = 0.0f;
                }
            }
        }
    }
    return out;
}

// Curve edges and border edges of the quad mesh.
BitArray lineEdges(const RawMesh& raw)
{
    const mesh::Adjacency& adjacency = raw.mesh.adjacency();
    BitArray lines = raw.curveEdges;
    const auto edgeCount = static_cast<std::size_t>(adjacency.edgeCount());
    if (lines.size() != edgeCount) {
        lines.resize(edgeCount);
    }
    for (const int e : adjacency.boundaryEdges()) {
        lines.set(static_cast<std::size_t>(e));
    }
    return lines;
}

// Clamps a length into [min, max] of edge e: below min gives min, else the
// smaller of max and the length (a length without number stays).
float clampLength(float d, float min, float max)
{
    if (min > d) {
        return min;
    }
    return max < d ? max : d;
}

bool isZero(Vec3f v) { return v.x == 0.0f && v.y == 0.0f && v.z == 0.0f; }

// z * f + p * (1 - f) per coordinate.
Vec3f pulled(Vec3f p, Vec3f z, float f)
{
    const float keep = 1.0f - f;
    return {z.x * f + p.x * keep, z.y * f + p.y * keep, z.z * f + p.z * keep};
}

// The nearer of the nearest points of z on the segments a-p and b-p (b on a tie).
Vec3f nearerOnSegments(Vec3f z, Vec3f a, Vec3f b, Vec3f p)
{
    const mesh::SegmentClosest onA = mesh::closestOnSegment(z, a, p);
    const mesh::SegmentClosest onB = mesh::closestOnSegment(z, b, p);
    return std::sqrt(onA.distanceSq) < std::sqrt(onB.distanceSq) ? onA.point : onB.point;
}

} // namespace

std::vector<float> FinalRelaxation::dihedralAngles(const mesh::Mesh& mesh)
{
    const mesh::Topology& topology = mesh.topology();
    const std::vector<Vec3f>& at = mesh.points();
    std::vector<Vec3f> normals(static_cast<std::size_t>(topology.faceCount()));
    std::vector<Vec3f> polygon;
    for (int f = 0; f < topology.faceCount(); ++f) {
        polygon.clear();
        for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
            polygon.push_back(at[static_cast<std::size_t>(topology.vertex(c))]);
        }
        const Vec3d n = prep::PolygonTriangulator::newellSum(polygon);
        const float length = std::sqrt(static_cast<float>((n.x * n.x + n.y * n.y) + n.z * n.z));
        const double l = static_cast<double>(length);
        if (1e-26 < l) {
            normals[static_cast<std::size_t>(f)] = {static_cast<float>(n.x / l), static_cast<float>(n.y / l),
                                                    static_cast<float>(n.z / l)};
        } else {
            normals[static_cast<std::size_t>(f)] = {static_cast<float>(n.x), static_cast<float>(n.y),
                                                    static_cast<float>(n.z)};
        }
    }
    return mesh::dihedralAngles(mesh.adjacency(), normals);
}

std::vector<Vec3f> FinalRelaxation::step(const mesh::Mesh& mesh, const StepInput& input)
{
    const mesh::Adjacency& adjacency = mesh.adjacency();
    const std::vector<Vec3f>& at = mesh.points();
    const std::vector<float> dihedral = dihedralAngles(mesh);
    const EdgeTargets& edges = input.edges;
    std::vector<Vec3f> out(at.size());
    auto onLine = [&](int e) { return input.lines.testBounded(static_cast<std::size_t>(e)); };
    // The mixes of the last point with four edges (they carry over where they are not computed).
    Vec3f aniso{};
    Vec3f iso{};
    for (int v = 0; v < mesh.pointCount(); ++v) {
        const auto i = static_cast<std::size_t>(v);
        const Vec3f p = at[i];
        out[i] = p;
        const float m = input.mobility[i];
        if (m == 0.0f) {
            continue;
        }
        float w = input.anisotropy[i];
        const auto pointEdges = adjacency.pointEdges(v);
        const std::size_t valence = pointEdges.size();
        if (valence == 0) {
            continue;
        }
        // Lengths, targets along the edge directions.
        auto length = [&](int e) {
            const auto k = static_cast<std::size_t>(e);
            return clampLength(mesh::length(at[static_cast<std::size_t>(adjacency.otherEnd(e, v))] - p),
                               edges.lengthMin[k], edges.lengthMax[k]);
        };
        auto target = [&](int e, float l) {
            const Vec3f r = edges.direction[static_cast<std::size_t>(e)];
            const Vec3f n = at[static_cast<std::size_t>(adjacency.otherEnd(e, v))];
            if (v == adjacency.edgeVertices(e)[0]) {
                return Vec3f{n.x - r.x * l, n.y - r.y * l, n.z - r.z * l};
            }
            return Vec3f{r.x * l + n.x, r.y * l + n.y, r.z * l + n.z};
        };
        auto pairTarget = [&](int ea, int eb) {
            const Vec3f qa = target(ea, length(ea));
            const Vec3f qb = target(eb, length(eb));
            return (qa + qb) * 0.5f;
        };
        auto neighbour = [&](int e) { return at[static_cast<std::size_t>(adjacency.otherEnd(e, v))]; };
        // The ordered ring around v when it has `cornerCount` corners and `edgeCount` edges.
        auto ring = [&](std::size_t cornerCount, std::size_t edgeCount) -> std::optional<mesh::VertexFan> {
            mesh::VertexFan fan = mesh::orderedVertexFan(adjacency, v);
            if (fan.corners.size() != adjacency.pointCorners(v).size() || fan.corners.size() != cornerCount ||
                fan.edges.size() != edgeCount) {
                return std::nullopt;
            }
            return fan;
        };

        if (valence == 4) {
            if (const std::optional<mesh::VertexFan> fan = ring(4, 5)) {
                const std::array<int, 4> e{fan->edges[0], fan->edges[1], fan->edges[2], fan->edges[3]};
                const bool b0 = onLine(e[0]);
                const bool b1 = onLine(e[1]);
                const bool b2 = onLine(e[2]);
                const bool b3 = onLine(e[3]);
                const int c02 = static_cast<int>(b0) + static_cast<int>(b2);
                const int c13 = static_cast<int>(b1) + static_cast<int>(b3);
                if (c02 == 1 || c13 == 1 || (c02 >= 1 && c13 >= 1)) {
                    continue;
                }
                bool useCross = true;
                if (c02 == 2) {
                    w = 1.0f;
                    useCross = !(mesh::wide(dihedral[static_cast<std::size_t>(e[0])] +
                                             dihedral[static_cast<std::size_t>(e[2])]) > creaseLimit);
                }
                if (c13 == 2) {
                    w = 1.0f;
                    useCross = !(mesh::wide(dihedral[static_cast<std::size_t>(e[1])] +
                                            dihedral[static_cast<std::size_t>(e[3])]) > creaseLimit);
                }
                if (w > 0.0f) {
                    const bool use02 = c13 == 2 ? useCross : true;
                    const bool use13 = c02 == 2 ? useCross : true;
                    auto usable = [&](int ea, int eb) {
                        return !isZero(edges.direction[static_cast<std::size_t>(ea)]) &&
                               !isZero(edges.direction[static_cast<std::size_t>(eb)]);
                    };
                    const Vec3f a = use02 && usable(e[0], e[2]) ? pairTarget(e[0], e[2]) : p;
                    const Vec3f b = use13 && usable(e[1], e[3]) ? pairTarget(e[1], e[3]) : p;
                    aniso = (a + b) * 0.5f;
                }
                if (w < 1.0f) {
                    auto between = [&](int ea, int eb) {
                        const float la = length(ea);
                        const float lb = length(eb);
                        const float f = la / (lb + la);
                        const Vec3f na = neighbour(ea);
                        const Vec3f nb = neighbour(eb);
                        const float g = 1.0f - f;
                        return Vec3f{na.x * g + nb.x * f, na.y * g + nb.y * f, na.z * g + nb.z * f};
                    };
                    iso = (between(e[0], e[2]) + between(e[1], e[3])) * 0.5f;
                }
                const float keep = 1.0f - w;
                Vec3f z{iso.x * keep + aniso.x * w, iso.y * keep + aniso.y * w, iso.z * keep + aniso.z * w};
                if (b0 && b2) {
                    z = nearerOnSegments(z, neighbour(e[0]), neighbour(e[2]), p);
                }
                if (b1 && b3) {
                    z = nearerOnSegments(z, neighbour(e[1]), neighbour(e[3]), p);
                }
                out[i] = pulled(p, z, m * pullFactor);
                continue;
            }
        } else if (valence == 3) {
            bool border = false;
            for (const int e : pointEdges) {
                border = border || adjacency.isBoundaryEdge(e);
            }
            if (border && w > 0.0f) {
                if (const std::optional<mesh::VertexFan> fan = ring(2, 3)) {
                    const int e0 = fan->edges[0];
                    const int e1 = fan->edges[1];
                    const int e2 = fan->edges[2];
                    const Vec3f a = isZero(edges.direction[static_cast<std::size_t>(e0)]) ||
                                            isZero(edges.direction[static_cast<std::size_t>(e2)])
                                        ? p
                                        : pairTarget(e0, e2);
                    if (onLine(e1)) {
                        continue;
                    }
                    // The inner edge with its plain length.
                    const float d1 = mesh::length(neighbour(e1) - p);
                    const Vec3f q1 = target(e1, d1);
                    Vec3f z = (a + q1) * 0.5f;
                    z = nearerOnSegments(z, neighbour(e0), neighbour(e2), p);
                    out[i] = pulled(p, z, m * pullFactor);
                    continue;
                }
            }
        }
        // Mean of the neighbours, along a line when the point lies on one.
        const auto neighbours = adjacency.pointNeighbours(v);
        Vec3f sum{};
        int lineCount = 0;
        std::array<Vec3f, 2> line{};
        for (std::size_t k = 0; k < valence; ++k) {
            const Vec3f q = at[static_cast<std::size_t>(neighbours[k])];
            if (onLine(pointEdges[k])) {
                if (lineCount < 2) {
                    line[static_cast<std::size_t>(lineCount)] = q;
                }
                ++lineCount;
            }
            sum = sum + q;
        }
        const float n = static_cast<float>(valence);
        Vec3f z{sum.x / n, sum.y / n, sum.z / n};
        if (lineCount == 2) {
            z = nearerOnSegments(z, line[0], line[1], p);
        } else if (lineCount != 0) {
            continue;
        }
        out[i] = pulled(p, z, m * meanFactor);
    }
    return out;
}

FinalRelaxation::EdgeTargets FinalRelaxation::edgeTargets(parts::Part& part)
{
    return EdgeTargetBuilder(part, *part.raw).build();
}

void FinalRelaxation::run(parts::Part& part, bool symmetry)
{
    if (!part.raw || !part.trace || !part.sizeField || !part.cutGraph || !part.crossField) {
        return;
    }
    RawMesh& raw = *part.raw;
    StepInput input;
    input.edges = edgeTargets(part);
    input.anisotropy = pointAnisotropy(part, raw);
    const BitArray fixed = fixedEdges(part, symmetry);
    input.mobility = pointMobility(part, raw, fixed);
    input.lines = lineEdges(raw);
    SurfaceWalker walker(part.mesh, fixed);
    std::vector<PointSource>& sources = raw.points.source;
    for (int round = 0; round < rounds; ++round) {
        std::vector<Vec3f> moved = step(raw.mesh, input);
        const std::vector<Vec3f>& old = raw.mesh.points();
        for (std::size_t i = 0; i < moved.size(); ++i) {
            const Vec3f shift = moved[i] - old[i];
            moved[i] = old[i];
            if (i < sources.size()) {
                walker.walk(sources[i], moved[i], shift);
            }
        }
        raw.mesh.swapPoints(moved);
    }
}

} // namespace areko::remesher::extract
