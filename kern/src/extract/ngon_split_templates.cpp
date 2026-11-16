#include "ngon_split_templates.h"

#include "prep/polygon_triangulator.h"

#include <bit>
#include <cstdint>
#include <utility>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

constexpr float offsetFactor = std::bit_cast<float>(0x3f19999au); // 0.6

} // namespace

void SplitTemplates::walk(const NgonPolygonList& list, int p, int from, int to, std::vector<int>& out)
{
    out.push_back(list.vertex(from));
    for (int i = list.next(p, from); i != to; i = list.next(p, i)) {
        out.push_back(list.vertex(i));
    }
    out.push_back(list.vertex(to));
}

void SplitTemplates::replace(NgonSplitWork& work, int p, const std::vector<std::vector<int>>& parts)
{
    work.corners.erase(work.list.begin(p), work.list.end(p));
    work.list.remove(p);
    for (const auto& part : parts) {
        work.list.append(part);
    }
}

void SplitTemplates::cut(NgonSplitWork& work, int p, int a, int b)
{
    if (a > b) {
        std::swap(a, b);
    }
    const auto q = work.list.polygon(p);
    std::vector<std::vector<int>> parts(2);
    parts[0].assign(q.begin(), q.begin() + a + 1);
    parts[0].insert(parts[0].end(), q.begin() + b, q.end());
    parts[1].assign(q.begin() + a, q.begin() + b + 1);
    replace(work, p, parts);
}

void SplitTemplates::singleCorner(NgonSplitWork& work, int p, int a, int b, int c)
{
    const NgonPolygonList& list = work.list;
    NgonNewPoints& points = work.points;
    const Vec3f pa = points.position(list.vertex(a));
    const Vec3f pb = points.position(list.vertex(b));
    const Vec3f pc = points.position(list.vertex(c));
    const int x = points.add(((pa - pb) + pb) + (pc - pb));
    points.addParent(list.vertex(b), -1.0f);
    points.addParent(list.vertex(a), 1.0f);
    points.addParent(list.vertex(c), 1.0f);
    points.endParents();
    const int m = list.wrap(p, b + list.size(p) / 2);
    std::vector<std::vector<int>> parts(3);
    walk(list, p, a, c, parts[0]);
    walk(list, p, c, m, parts[1]);
    walk(list, p, m, a, parts[2]);
    for (auto& part : parts) {
        part.push_back(x);
    }
    replace(work, p, parts);
}

void SplitTemplates::threeCorners(NgonSplitWork& work, int p, int v0, int v1, int v2, int v3, int v4)
{
    const NgonPolygonList& list = work.list;
    NgonNewPoints& points = work.points;
    auto at = [&](int i) { return points.position(list.vertex(i)); };
    const float l1 = length(at(v0) - at(v1));
    const float l2 = length(at(v4) - at(v3));
    Vec3f e1 = at(v3) - at(v2);
    Vec3f e2 = at(v1) - at(v2);
    normalizeInPlace(e1);
    normalizeInPlace(e2);
    const int x = points.add((at(v2) + e1 * l1) + e2 * l2);
    points.addParent(list.vertex(v2), 1.0f);
    points.addParent(list.vertex(v0), 1.0f);
    points.addParent(list.vertex(v1), -1.0f);
    points.addParent(list.vertex(v4), 1.0f);
    points.addParent(list.vertex(v3), -1.0f);
    points.endParents();
    std::vector<std::vector<int>> parts(3);
    parts[0] = {list.vertex(v0), list.vertex(v1), list.vertex(v2), x};
    parts[1] = {list.vertex(v2), list.vertex(v3), list.vertex(v4), x};
    walk(list, p, v4, v0, parts[2]);
    parts[2].push_back(x);
    replace(work, p, parts);
}

void SplitTemplates::cornerStraightsCorner(NgonSplitWork& work, int p, int k, int m)
{
    const NgonPolygonList& list = work.list;
    NgonNewPoints& points = work.points;
    auto at = [&](int i) { return points.position(list.vertex(i)); };
    const int kp = list.prev(p, k);
    const int kn = list.next(p, k);
    const int mp = list.prev(p, m);
    const int mn = list.next(p, m);
    const int x1 = points.add((at(kp) - at(k)) * offsetFactor + at(kn));
    points.addParent(list.vertex(kn), 1.0f);
    points.addParent(list.vertex(kp), offsetFactor);
    points.addParent(list.vertex(k), -offsetFactor);
    points.endParents();
    const int x2 = points.add((at(mn) - at(m)) * offsetFactor + at(mp));
    points.addParent(list.vertex(mp), 1.0f);
    points.addParent(list.vertex(mn), offsetFactor);
    points.addParent(list.vertex(m), -offsetFactor);
    points.endParents();
    std::vector<std::vector<int>> parts(4);
    walk(list, p, kn, mp, parts[0]);
    parts[0].push_back(x2);
    parts[0].push_back(x1);
    parts[1] = {list.vertex(m), list.vertex(mn), x2, list.vertex(mp)};
    parts[2] = {x1, x2};
    walk(list, p, mn, kp, parts[2]);
    parts[3] = {list.vertex(k), list.vertex(kn), x1, list.vertex(kp)};
    replace(work, p, parts);
}

void SplitTemplates::fan(NgonSplitWork& work, int p, int a, int b, int c)
{
    const auto q = work.list.polygon(p);
    const int n = static_cast<int>(q.size());
    const int s = work.list.begin(p);
    a -= s;
    b -= s;
    c -= s;
    const int db = b <= a ? b - a + n : b - a;
    const int dc = c <= a ? c - a + n : c - a;
    if (dc < db) {
        std::swap(b, c);
    }
    // Local corners from `from` to `to` (both included), round the end if needed.
    auto range = [&](int from, int to, std::vector<int>& out) {
        const int f = (to + 1) % n;
        if (f <= from) {
            out.insert(out.end(), q.begin() + from, q.end());
            out.insert(out.end(), q.begin(), q.begin() + f);
        } else {
            out.insert(out.end(), q.begin() + from, q.begin() + to + 1);
        }
    };
    std::vector<std::vector<int>> parts(3);
    range(a, b, parts[0]);
    range(b, c, parts[1]);
    parts[1].push_back(q[static_cast<std::size_t>(a)]);
    range(c, a, parts[2]);
    replace(work, p, parts);
}

void SplitTemplates::triangulate(NgonSplitWork& work, int p)
{
    const auto q = work.list.polygon(p);
    const int n = static_cast<int>(q.size());
    std::vector<Vec3f> polygon;
    polygon.reserve(q.size());
    for (const int v : q) {
        polygon.push_back(work.points.position(v));
    }
    prep::LocalTriangles triangles;
    if (!prep::PolygonTriangulator::tessellate(polygon, triangles, true)) {
        triangles.clear();
        prep::PolygonTriangulator::fan(n, 0, triangles);
    }
    std::vector<std::vector<int>> parts;
    parts.reserve(triangles.size());
    for (const auto& t : triangles) {
        parts.push_back({q[static_cast<std::size_t>(t[0])], q[static_cast<std::size_t>(t[1])],
                         q[static_cast<std::size_t>(t[2])]});
    }
    replace(work, p, parts);
}

} // namespace areko::remesher::extract
