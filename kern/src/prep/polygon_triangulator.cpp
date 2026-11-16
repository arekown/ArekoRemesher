#include "polygon_triangulator.h"

#include <CDT.h>
#include <GL/glu.h>

#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <utility>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

// Axis e for the plane basis u = e x n: the axis along the smallest normal component.
Vec3d basisAxis(const Vec3d& n)
{
    const double ax = std::fabs(n.x);
    const double ay = std::fabs(n.y);
    const double az = std::fabs(n.z);
    if (ax <= ay && ax <= az) {
        return {1.0, 0.0, 0.0};
    }
    if (ay < ax && ay < az) {
        return {0.0, 1.0, 0.0};
    }
    return {0.0, 0.0, 1.0};
}

// Plane axis u = e x n of the tessellator, not normalized, written out per
// axis: e = x gives (0, -n.z, n.y), e = y (n.z, 0, -n.x), e = z (-n.y, n.x, 0).
// e is the axis along the smallest normal component; when |n.x| or |n.y| is
// NaN, e = x unless |n.x| > |n.z| (then z), and e = y needs both comparisons
// to hold.
Vec3f tessellatorAxis(const Vec3f& n)
{
    const float ax = std::fabs(n.x);
    const float ay = std::fabs(n.y);
    const float az = std::fabs(n.z);
    if (!(ax > ay) && !(ax > az)) {
        return {0.0f, -n.z, n.y};
    }
    if (ay < ax && ay < az) {
        return {n.z, 0.0f, -n.x};
    }
    return {-n.y, n.x, 0.0f};
}

// Newell normal: sums over the edges (i, i+1), the closing edge last.
template <class V>
V newell(std::span<const Vec3f> polygon)
{
    using T = decltype(V{}.x);
    V n{};
    const std::size_t count = polygon.size();
    for (std::size_t i = 0; i < count; ++i) {
        const Vec3f& a = polygon[i];
        const Vec3f& b = polygon[(i + 1) % count];
        n.x += (T(a.y) - T(b.y)) * (T(a.z) + T(b.z));
        n.y += (T(a.z) - T(b.z)) * (T(a.x) + T(b.x));
        n.z += (T(a.x) - T(b.x)) * (T(a.y) + T(b.y));
    }
    return n;
}

using Point2 = std::pair<double, double>;

// Moves every point that equals an earlier one along the sum of its two edge
// vectors, by the first factor 1e-15, 1e-14, ... 0.1 that changes it.
// False when equal points remain.
bool separateDuplicates(std::vector<Point2>& points)
{
    const std::size_t n = points.size();
    std::vector<int> duplicateOf(n, -1);
    std::map<Point2, int> first;
    for (std::size_t i = 0; i < n; ++i) {
        const auto [it, inserted] = first.emplace(points[i], static_cast<int>(i));
        if (!inserted) {
            duplicateOf[i] = it->second;
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (duplicateOf[i] < 0) {
            continue;
        }
        const Point2 p = points[i];
        const Point2 next = points[(i + 1) % n];
        const Point2 prev = points[(i + n - 1) % n];
        const double dx = (next.first - p.first) + (prev.first - p.first);
        const double dy = (next.second - p.second) + (prev.second - p.second);
        double k = 1e-15;
        for (int step = 0; step < 15; ++step, k *= 10.0) {
            const Point2 moved{dx * k + p.first, dy * k + p.second};
            if (moved.first != p.first || moved.second != p.second) {
                points[i] = moved;
                break;
            }
        }
    }
    std::set<Point2> distinct(points.begin(), points.end());
    return distinct.size() == n;
}

// Tessellator callbacks collect the vertex numbers of every primitive.
struct TessState {
    LocalTriangles triangles;
    std::vector<int> primitive;
    GLenum type = 0;
    bool combined = false;
    bool failed = false;
};

void GLAPIENTRY onBegin(GLenum type, void* data)
{
    auto& state = *static_cast<TessState*>(data);
    state.type = type;
    state.primitive.clear();
}

void GLAPIENTRY onVertex(void* vertex, void* data)
{
    static_cast<TessState*>(data)->primitive.push_back(static_cast<int>(reinterpret_cast<std::intptr_t>(vertex)));
}

void GLAPIENTRY onEnd(void* data)
{
    auto& state = *static_cast<TessState*>(data);
    const auto& v = state.primitive;
    const std::size_t m = v.size();
    switch (state.type) {
    case GL_TRIANGLES:
        for (std::size_t j = 0; j + 2 < m; j += 3) {
            state.triangles.push_back({v[j], v[j + 1], v[j + 2]});
        }
        break;
    case GL_TRIANGLE_STRIP:
        for (std::size_t j = 0; j + 2 < m; ++j) {
            if (j % 2 == 0) {
                state.triangles.push_back({v[j], v[j + 1], v[j + 2]});
            } else {
                state.triangles.push_back({v[j], v[j + 2], v[j + 1]});
            }
        }
        break;
    case GL_TRIANGLE_FAN:
        for (std::size_t j = 1; j + 1 < m; ++j) {
            state.triangles.push_back({v[0], v[j], v[j + 1]});
        }
        break;
    default:
        break;
    }
}

void GLAPIENTRY onError(GLenum, void* data)
{
    static_cast<TessState*>(data)->failed = true;
}

void GLAPIENTRY onCombine(GLdouble*, void**, GLfloat*, void** out, void* data)
{
    static_cast<TessState*>(data)->combined = true;
    *out = reinterpret_cast<void*>(std::intptr_t{-1});
}

// Owns one tessellator object.
class Tessellator {
public:
    Tessellator() : tess_(gluNewTess()) {}
    ~Tessellator()
    {
        if (tess_ != nullptr) {
            gluDeleteTess(tess_);
        }
    }
    Tessellator(const Tessellator&) = delete;
    Tessellator& operator=(const Tessellator&) = delete;
    [[nodiscard]] GLUtesselator* get() const noexcept { return tess_; }

private:
    GLUtesselator* tess_;
};

} // namespace

void PolygonTriangulator::fan(int n, int start, LocalTriangles& out)
{
    for (int k = 0; k < n; ++k) {
        const int k1 = (k + 1) % n;
        if (k != start && k1 != start) {
            out.push_back({start, k, k1});
        }
    }
}

bool PolygonTriangulator::delaunay(std::span<const Vec3f> polygon, LocalTriangles& out)
{
    const std::size_t n = polygon.size();
    Vec3d normal = newell<Vec3d>(polygon);
    normalizeInPlace(normal);
    const Vec3d u = cross(basisAxis(normal), normal);
    const Vec3d w = cross(normal, u);
    std::vector<Point2> points(n);
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3d p{polygon[i].x, polygon[i].y, polygon[i].z};
        points[i] = {(u.x * p.x + u.y * p.y) + u.z * p.z, (p.x * w.x + p.y * w.y) + p.z * w.z};
    }
    if (!separateDuplicates(points)) {
        return false;
    }
    try {
        CDT::Triangulation<double> cdt(CDT::VertexInsertionOrder::Auto, CDT::IntersectingConstraintEdges::Resolve, 0.0);
        std::vector<CDT::V2d<double>> vertices;
        vertices.reserve(n);
        for (const auto& [x, y] : points) {
            vertices.push_back(CDT::V2d<double>::make(x, y));
        }
        cdt.insertVertices(vertices);
        std::vector<CDT::Edge> edges;
        edges.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            edges.emplace_back(static_cast<CDT::VertInd>(i), static_cast<CDT::VertInd>((i + 1) % n));
        }
        cdt.insertEdges(edges);
        cdt.eraseOuterTriangles();
        if (cdt.triangles.size() != n - 2 || cdt.vertices.size() > n) {
            return false;
        }
        for (const CDT::Triangle& t : cdt.triangles) {
            out.push_back({static_cast<int>(t.vertices[0]), static_cast<int>(t.vertices[1]), static_cast<int>(t.vertices[2])});
        }
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

Vec3d PolygonTriangulator::newellSum(std::span<const Vec3f> polygon)
{
    return newell<Vec3d>(polygon);
}

Vec3f PolygonTriangulator::planeNormal(std::span<const Vec3f> polygon)
{
    Vec3f normal = newell<Vec3f>(polygon);
    normalizeInPlace(normal);
    return normal;
}

bool PolygonTriangulator::tessellate(std::span<const Vec3f> polygon, LocalTriangles& out, bool rejectOnError)
{
    const std::size_t n = polygon.size();
    const Vec3f normal = planeNormal(polygon);
    const Vec3f u = tessellatorAxis(normal);
    const Vec3f w = cross(normal, u);
    std::vector<std::array<GLdouble, 3>> coords(n);
    for (std::size_t i = 0; i < n; ++i) {
        coords[i] = {static_cast<double>(dot(u, polygon[i])), static_cast<double>(dot(w, polygon[i])), 0.0};
    }
    Tessellator tess;
    if (tess.get() == nullptr) {
        return false;
    }
    TessState state;
    gluTessProperty(tess.get(), GLU_TESS_WINDING_RULE, GLU_TESS_WINDING_NONZERO);
    gluTessProperty(tess.get(), GLU_TESS_BOUNDARY_ONLY, GL_FALSE);
    gluTessCallback(tess.get(), GLU_TESS_BEGIN_DATA, reinterpret_cast<_GLUfuncptr>(&onBegin));
    gluTessCallback(tess.get(), GLU_TESS_VERTEX_DATA, reinterpret_cast<_GLUfuncptr>(&onVertex));
    gluTessCallback(tess.get(), GLU_TESS_END_DATA, reinterpret_cast<_GLUfuncptr>(&onEnd));
    gluTessCallback(tess.get(), GLU_TESS_ERROR_DATA, reinterpret_cast<_GLUfuncptr>(&onError));
    gluTessCallback(tess.get(), GLU_TESS_COMBINE_DATA, reinterpret_cast<_GLUfuncptr>(&onCombine));
    gluTessBeginPolygon(tess.get(), &state);
    gluTessBeginContour(tess.get());
    for (std::size_t i = 0; i < n; ++i) {
        gluTessVertex(tess.get(), coords[i].data(), reinterpret_cast<void*>(static_cast<std::intptr_t>(i)));
    }
    gluTessEndContour(tess.get());
    gluTessEndPolygon(tess.get());
    if (state.combined || (rejectOnError && state.failed) || state.triangles.size() != n - 2) {
        return false;
    }
    out = std::move(state.triangles);
    return true;
}

LocalTriangles PolygonTriangulator::triangulate(std::span<const Vec3f> polygon)
{
    LocalTriangles out;
    const int n = static_cast<int>(polygon.size());
    if (n < 3) {
        return out;
    }
    if (delaunay(polygon, out)) {
        return out;
    }
    out.clear();
    if (tessellate(polygon, out)) {
        return out;
    }
    out.clear();
    fan(n, 0, out);
    return out;
}

} // namespace areko::remesher::prep
