#include "normals.h"

#include "attribute_maps.h"
#include "bit_array.h"
#include "mesh.h"
#include "vertex_fan.h"

#include <bit>
#include <cmath>

namespace areko::remesher::mesh {

void floatTriangleNormal(Vec3f a, Vec3f b, Vec3f c, Vec3f& normal, float& doubleArea)
{
    const TriangleNormal n = triangleNormal(a, b, c);
    normal = {static_cast<float>(n.unit.x), static_cast<float>(n.unit.y), static_cast<float>(n.unit.z)};
    doubleArea = static_cast<float>(n.length);
}

std::vector<Vec3f> triangleCentres(const Mesh& mesh)
{
    const Topology& topology = mesh.topology();
    const std::vector<Vec3f>& points = mesh.points();
    std::vector<Vec3f> out(static_cast<std::size_t>(mesh.faceCount()));
    for (std::size_t t = 0; t < out.size(); ++t) {
        const int c = 3 * static_cast<int>(t);
        const Vec3d p0 = wide(points[static_cast<std::size_t>(topology.vertex(c))]);
        const Vec3d p1 = wide(points[static_cast<std::size_t>(topology.vertex(c + 1))]);
        const Vec3d p2 = wide(points[static_cast<std::size_t>(topology.vertex(c + 2))]);
        out[t] = {static_cast<float>(((p0.x + p1.x) + p2.x) / 3.0), static_cast<float>(((p0.y + p1.y) + p2.y) / 3.0),
                  static_cast<float>(((p0.z + p1.z) + p2.z) / 3.0)};
    }
    return out;
}

namespace {

// Triangle over three corners; zero normal and area when a corner has no point.
void cornerTriangle(const Topology& topo, const std::vector<Vec3f>& points, int c0, int c1, int c2, Vec3f& normal,
                    float& doubleArea)
{
    normal = {};
    doubleArea = 0.0f;
    const int v0 = topo.vertex(c0), v1 = topo.vertex(c1), v2 = topo.vertex(c2);
    if (v0 >= 0 && v1 >= 0 && v2 >= 0) {
        floatTriangleNormal(points[static_cast<std::size_t>(v0)], points[static_cast<std::size_t>(v1)],
                            points[static_cast<std::size_t>(v2)], normal, doubleArea);
    }
}

} // namespace

TriangleGeometry::TriangleGeometry(const Mesh& mesh)
{
    const Topology& topo = mesh.topology();
    const auto& points = mesh.points();
    float sum = 0.0f;
    auto add = [&](int c0, int c1, int c2) {
        Vec3f n;
        float a2 = 0.0f;
        cornerTriangle(topo, points, c0, c1, c2, n, a2);
        normals_.push_back(n);
        doubleAreas_.push_back(a2);
        sum += a2;
    };
    const int faces = topo.faceCount();
    for (int f = 0; f < faces; ++f) {
        const int begin = topo.faceBegin(f);
        const int end = topo.faceEnd(f);
        for (int c = begin + 1; c + 1 < end; ++c) {
            add(begin, c, c + 1);
        }
    }
    totalArea_ = 0.5f * sum;
}

float TriangleGeometry::smallAreaLimit() const noexcept
{
    constexpr float minimum = std::bit_cast<float>(0x2b8cbcccu); // 1e-12
    const auto limit = static_cast<float>(static_cast<double>(totalArea_) * 1e-8);
    return static_cast<double>(limit) < 1e-12 ? minimum : limit;
}

void TriangleGeometry::replaceDegenerateNormals(const Mesh& mesh, const Adjacency& adjacency, int mode)
{
    if (mesh.topology().mode() != FaceMode::Triangles) {
        return;
    }
    const std::size_t count = normals_.size();
    const float limit = smallAreaLimit();
    BitArray degenerate(count);
    for (std::size_t t = 0; t < count; ++t) {
        if (!(doubleAreas_[t] > limit)) {
            degenerate.set(t);
        }
    }
    int rest = static_cast<int>(degenerate.count());
    if (rest == 0 || mode != 2) {
        return;
    }
    const auto& points = mesh.points();
    auto edgeLength = [&](int e) {
        const auto ends = adjacency.edgeVertices(e);
        return length(points[static_cast<std::size_t>(ends[0])] - points[static_cast<std::size_t>(ends[1])]);
    };

    for (int rounds = 20;;) {
        const int start = rest;
        for (std::ptrdiff_t ti = degenerate.findNext(0); ti >= 0; ti = degenerate.findNext(static_cast<std::size_t>(ti) + 1)) {
            const int t = static_cast<int>(ti);
            int validNb[3];
            float degenerateLen[3], validLen[3];
            int d = 0, g = 0;
            for (int i = 0; i < 3; ++i) {
                const int c = 3 * t + i;
                const int e = adjacency.cornerEdge(c);
                if (e < 0) {
                    continue;
                }
                const auto corners = adjacency.edgeCorners(e);
                if (corners.size() != 2) {
                    continue;
                }
                int other;
                if (c == corners[0]) {
                    other = corners[1];
                } else if (c == corners[1]) {
                    other = corners[0];
                } else {
                    continue;
                }
                if (other < 0) {
                    continue;
                }
                const int u = other / 3;
                const float len = edgeLength(e);
                if (static_cast<std::size_t>(u) < count && degenerate.test(static_cast<std::size_t>(u))) {
                    degenerateLen[d++] = len;
                } else {
                    validNb[g] = u;
                    validLen[g++] = len;
                }
            }
            int source = -1;
            bool replaced = false;
            if (d == 2 && g == 1) {
                if (degenerateLen[0] > validLen[0] || degenerateLen[1] > validLen[0]) {
                    source = validNb[0];
                }
            } else if (d == 1 && g == 2) {
                const float p = dot(normals_[static_cast<std::size_t>(validNb[0])], normals_[static_cast<std::size_t>(validNb[1])]);
                if (!(0.75f > p)) {
                    if (degenerateLen[0] > validLen[0] && degenerateLen[0] > validLen[1]) {
                        source = validLen[0] > validLen[1] ? validNb[0] : validNb[1];
                    } else {
                        source = validLen[0] > validLen[1] ? validNb[1] : validNb[0];
                    }
                }
            } else if (g == 3) {
                int order[3] = {0, 1, 2};
                for (int i = 1; i < 3; ++i) { // stable insertion sort, strict comparison
                    const int key = order[i];
                    int j = i - 1;
                    while (j >= 0 && validLen[key] < validLen[order[j]]) {
                        order[j + 1] = order[j];
                        --j;
                    }
                    order[j + 1] = key;
                }
                const Vec3f a = normals_[static_cast<std::size_t>(validNb[order[0]])];
                const Vec3f b = normals_[static_cast<std::size_t>(validNb[order[1]])];
                if (!(0.75f > dot(a, b))) {
                    Vec3f v = a + b;
                    normalizeInPlace(v);
                    normals_[static_cast<std::size_t>(t)] = v;
                    replaced = true;
                }
            } else if (g >= 1) {
                source = validNb[0];
            }
            if (source >= 0) {
                normals_[static_cast<std::size_t>(t)] = normals_[static_cast<std::size_t>(source)];
                replaced = true;
            }
            if (replaced) {
                degenerate.set(static_cast<std::size_t>(t), false);
                --rest;
            }
        }
        if (rest == start || rest <= 0 || --rounds == 0) {
            return;
        }
    }
}

void TriangleGeometry::updateAt(const Mesh& mesh, std::span<const int> points)
{
    const Topology& topo = mesh.topology();
    if (topo.mode() != FaceMode::Triangles || points.empty()) {
        return;
    }
    const Adjacency& adjacency = mesh.adjacency();
    const auto& positions = mesh.points();
    BitArray touched(normals_.size());
    for (const int v : points) {
        for (const int c : adjacency.pointCorners(v)) {
            touched.set(static_cast<std::size_t>(c / 3));
        }
    }
    for (std::ptrdiff_t ti = touched.findNext(0); ti >= 0; ti = touched.findNext(static_cast<std::size_t>(ti) + 1)) {
        const auto t = static_cast<std::size_t>(ti);
        const int c = 3 * static_cast<int>(ti);
        cornerTriangle(topo, positions, c, c + 1, c + 2, normals_[t], doubleAreas_[t]);
    }
    float sum = 0.0f;
    for (const float a2 : doubleAreas_) {
        sum += a2;
    }
    totalArea_ = 0.5f * sum;
}

namespace {

// Smallest cosine between neighbouring triangles of one group (about 50 degrees).
constexpr float groupCosine = std::bit_cast<float>(0x3f248f05u);

} // namespace

SplitNormals::SplitNormals(const Mesh& mesh)
{
    const auto& given = mesh.customNormals();
    if (given && given->cornerNormal.size() == static_cast<std::size_t>(mesh.cornerCount())) {
        cornerNormal_ = given->cornerNormal;
        normals_ = given->normals;
        return;
    }
    if (mesh.topology().mode() == FaceMode::Triangles) {
        build(mesh);
    }
}

void SplitNormals::build(const Mesh& mesh)
{
    const Adjacency& adjacency = mesh.adjacency();
    const TriangleGeometry& triangles = mesh.triangleGeometry();
    const auto& faceNormals = triangles.normals();
    const auto& doubleAreas = triangles.doubleAreas();
    const auto* sharpMap = mesh.maps().findAs<BoolMap>(map_names::sharpEdges, ElementKind::HalfEdge);
    const BitArray* sharp = sharpMap != nullptr ? &sharpMap->bits() : nullptr;
    auto normal = [&](int c) { return faceNormals[static_cast<std::size_t>(adjacency.cornerFace(c))]; };
    auto weighted = [&](int c) {
        const auto t = static_cast<std::size_t>(adjacency.cornerFace(c));
        return faceNormals[t] * doubleAreas[t];
    };
    auto isSharp = [&](int c) { return sharp != nullptr && sharp->testBounded(static_cast<std::size_t>(c)); };
    auto slot = [&](int c) -> int& { return cornerNormal_[static_cast<std::size_t>(c)]; };

    const auto cornerCount = static_cast<std::size_t>(mesh.cornerCount());
    cornerNormal_.assign(cornerCount, -1);
    normals_.reserve(cornerCount / 2);
    for (int v = 0; v < mesh.pointCount(); ++v) {
        const std::vector<int> fan = orderedVertexFan(adjacency, v).corners;
        if (!fan.empty()) {
            slot(fan[0]) = static_cast<int>(normals_.size());
            normals_.push_back(weighted(fan[0]));
            for (std::size_t i = 1; i < fan.size(); ++i) {
                const int c = fan[i];
                // A NaN cosine keeps the group.
                if (isSharp(c) || groupCosine > dot(normal(fan[i - 1]), normal(c))) {
                    slot(c) = static_cast<int>(normals_.size());
                    normals_.push_back(weighted(c));
                } else {
                    normals_.back() = weighted(c) + normals_.back();
                    slot(c) = slot(fan[i - 1]);
                }
            }
            const int first = fan.front();
            const int last = fan.back();
            const float cosine = dot(normal(last), normal(first));
            // A NaN cosine joins with a sharp edge map and stays apart without one.
            const bool join = sharp != nullptr ? !(groupCosine > cosine) && !isSharp(first) : cosine >= groupCosine;
            const int lastGroup = slot(last);
            const int firstGroup = slot(first);
            if (join && lastGroup != firstGroup) {
                normals_.pop_back();
                for (std::size_t i = fan.size(); i-- > 0 && slot(fan[i]) == lastGroup;) {
                    slot(fan[i]) = firstGroup;
                }
            }
        }
        int placeholder = -1;
        for (const int c : adjacency.pointCorners(v)) {
            if (slot(c) < 0) {
                if (placeholder < 0) {
                    placeholder = static_cast<int>(normals_.size());
                    normals_.push_back({1.0f, 0.0f, 0.0f});
                }
                slot(c) = placeholder;
            }
        }
    }
    for (Vec3f& n : normals_) {
        normalizeInPlace(n);
    }
}

Vec3f SplitNormals::atCorner(int corner) const
{
    if (corner < 0 || static_cast<std::size_t>(corner) >= cornerNormal_.size()) {
        return {};
    }
    const int i = cornerNormal_[static_cast<std::size_t>(corner)];
    return i >= 0 && static_cast<std::size_t>(i) < normals_.size() ? normals_[static_cast<std::size_t>(i)] : Vec3f{};
}

Vec3f SplitNormals::edgeNormal(const Adjacency& adjacency, int e) const
{
    const auto corners = adjacency.edgeCorners(e);
    if (corners.size() > 2) {
        return {};
    }
    Vec3f sum{};
    for (const int c : corners) {
        sum = sum + atCorner(c);
        sum = sum + atCorner(adjacency.next(c));
    }
    normalizeInPlace(sum);
    return sum;
}

std::vector<float> dihedralAngles(const Adjacency& adjacency, const std::vector<Vec3f>& faceNormals)
{
    std::vector<float> angles(static_cast<std::size_t>(adjacency.edgeCount()), 0.0f);
    for (int e = 0; e < adjacency.edgeCount(); ++e) {
        const auto corners = adjacency.edgeCorners(e);
        if (corners.size() != 2) {
            continue;
        }
        const float p = dot(faceNormals[static_cast<std::size_t>(adjacency.cornerFace(corners[0]))],
                            faceNormals[static_cast<std::size_t>(adjacency.cornerFace(corners[1]))]);
        angles[static_cast<std::size_t>(e)] = static_cast<float>(angleOfCosine(static_cast<double>(p)));
    }
    return angles;
}

} // namespace areko::remesher::mesh
