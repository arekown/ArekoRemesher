#include "ngon_parity_balance.h"

#include "util/index_queue.h"
#include "util/introsort.h"

#include <algorithm>
#include <limits>
#include <map>
#include <utility>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

constexpr int largestEndFace = 7;
constexpr float innerShare = 0.4f;
constexpr float cornerShare = 0.2f;
constexpr float innerDivisor = 2.5f;

// First corner of edge e (edge corner order) that belongs to face f.
int cornerIn(const Adjacency& adjacency, int edge, int face)
{
    for (const int c : adjacency.edgeCorners(edge)) {
        if (adjacency.cornerFace(c) == face) {
            return c;
        }
    }
    return -1;
}

// Point ids of the corners from `from` round the face to `to` (both included).
void cornerRun(const Adjacency& adjacency, int from, int to, std::vector<int>& out)
{
    for (int c = from;; c = adjacency.next(c)) {
        out.push_back(adjacency.vertex(c));
        if (c == to) {
            break;
        }
    }
}

// New point in the middle of a and b, with both as parents.
int middlePoint(NgonNewPoints& points, int a, int b)
{
    const int id = points.add((points.position(a) + points.position(b)) * 0.5f);
    points.addParent(a, 0.5f);
    points.addParent(b, 0.5f);
    points.endParents();
    return id;
}

} // namespace

ParityBalance::Search ParityBalance::search(const Mesh& mesh, const BitArray& faces)
{
    const Topology& topology = mesh.topology();
    const Adjacency& adjacency = mesh.adjacency();
    const int faceCount = topology.faceCount();
    Search found;
    found.distance.assign(static_cast<std::size_t>(faceCount), std::numeric_limits<float>::max());
    found.source.assign(static_cast<std::size_t>(faceCount), -1);
    found.entry.assign(static_cast<std::size_t>(faceCount), -1);
    std::vector<int> odd;
    faces.forEachSet([&](std::size_t f) {
        if (topology.faceSize(static_cast<int>(f)) % 2 != 0) {
            odd.push_back(static_cast<int>(f));
        }
    });
    util::IndexQueue queue(static_cast<std::size_t>(faceCount));
    for (std::size_t i = 0; i < odd.size(); ++i) {
        const auto f = static_cast<std::size_t>(odd[i]);
        if (found.distance[f] == std::numeric_limits<float>::max()) {
            found.distance[f] = 0.0f;
            found.source[f] = static_cast<int>(i);
            found.entry[f] = -1;
            queue.pushBack(static_cast<int>(f));
        }
    }
    while (!queue.empty()) {
        const int f = queue.pop();
        for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
            const int e = adjacency.cornerEdge(c);
            if (e == -1) {
                continue;
            }
            const int g = adjacency.otherFace(e, f);
            if (g == -1 || !faces.testBounded(static_cast<std::size_t>(g))) {
                continue;
            }
            const float reach = found.distance[static_cast<std::size_t>(f)] + 1.0f;
            if (reach < found.distance[static_cast<std::size_t>(g)]) {
                found.distance[static_cast<std::size_t>(g)] = reach;
                found.source[static_cast<std::size_t>(g)] = found.source[static_cast<std::size_t>(f)];
                found.entry[static_cast<std::size_t>(g)] = e;
                queue.pushBack(g);
            }
        }
    }
    return found;
}

void ParityBalance::path(const Mesh& mesh, const Search& found, int start, std::vector<int>& faces,
                         std::vector<int>& edges)
{
    const Adjacency& adjacency = mesh.adjacency();
    faces.assign(1, start);
    edges.clear();
    int x = start;
    while (x >= 0 && found.entry[static_cast<std::size_t>(x)] != -1) {
        const int e = found.entry[static_cast<std::size_t>(x)];
        edges.push_back(e);
        x = adjacency.otherFace(e, x);
        faces.push_back(x);
    }
}

void ParityBalance::run(SplitRound& round, const BitArray& faces)
{
    const Mesh& mesh = round.mesh();
    if (mesh.cornerCount() == 0) {
        return;
    }
    const Topology& topology = mesh.topology();
    const Adjacency& adjacency = mesh.adjacency();
    round.clearDone();
    const Search found = search(mesh, faces);

    // Cheapest edge per pair of touching regions; on equal cost the first edge stays.
    std::map<std::pair<int, int>, int> pairIndex;
    std::vector<int> pairEdge;
    std::vector<float> pairCost;
    for (int e = 0; e < adjacency.edgeCount(); ++e) {
        const auto pair = adjacency.innerEdgeFaces(e);
        if (!pair) {
            continue;
        }
        const auto a = static_cast<std::size_t>((*pair)[0]);
        const auto b = static_cast<std::size_t>((*pair)[1]);
        const int qa = found.source[a];
        const int qb = found.source[b];
        if (qa == -1 || qb == -1 || qa == qb) {
            continue;
        }
        const float cost = (found.distance[a] + found.distance[b]) + 1.0f;
        const auto key = std::minmax(qa, qb);
        const auto [it, inserted] = pairIndex.emplace(key, static_cast<int>(pairEdge.size()));
        if (inserted) {
            pairEdge.push_back(e);
            pairCost.push_back(cost);
        } else if (cost < pairCost[static_cast<std::size_t>(it->second)]) {
            pairEdge[static_cast<std::size_t>(it->second)] = e;
            pairCost[static_cast<std::size_t>(it->second)] = cost;
        }
    }
    std::vector<int> order(pairEdge.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        order[i] = static_cast<int>(i);
    }
    util::Introsort::sort(order.begin(), order.end(), [&](int x, int y) {
        return pairCost[static_cast<std::size_t>(x)] < pairCost[static_cast<std::size_t>(y)];
    });

    std::vector<int> firstFaces;
    std::vector<int> firstEdges;
    std::vector<int> secondFaces;
    std::vector<int> secondEdges;
    std::vector<int> pathFaces;
    std::vector<int> pathEdges;
    for (const int k : order) {
        const int e = pairEdge[static_cast<std::size_t>(k)];
        const auto corners = adjacency.edgeCorners(e);
        path(mesh, found, adjacency.cornerFace(corners[0]), firstFaces, firstEdges);
        path(mesh, found, adjacency.cornerFace(corners[1]), secondFaces, secondEdges);
        pathFaces.assign(firstFaces.rbegin(), firstFaces.rend());
        pathFaces.insert(pathFaces.end(), secondFaces.begin(), secondFaces.end());
        pathEdges.assign(firstEdges.rbegin(), firstEdges.rend());
        pathEdges.push_back(e);
        pathEdges.insert(pathEdges.end(), secondEdges.begin(), secondEdges.end());
        if (topology.faceSize(pathFaces.back()) <= largestEndFace &&
            topology.faceSize(pathFaces.front()) <= largestEndFace) {
            splitPath(round, pathFaces, pathEdges);
        }
    }
}

void ParityBalance::splitPath(SplitRound& round, const std::vector<int>& faces, const std::vector<int>& edges)
{
    const std::size_t n = faces.size();
    if (n < 2) {
        return;
    }
    for (const int f : faces) {
        if (round.done(f)) {
            return;
        }
    }
    int middle = -1;
    endFace(round, edges[0], faces[0], middle);
    for (std::size_t i = 1; i + 1 < n; ++i) {
        const int entryMiddle = middle;
        middle = -1;
        innerFace(round, edges[i - 1], edges[i], faces[i], entryMiddle, middle);
    }
    endFace(round, edges[n - 2], faces[n - 1], middle);
}

void ParityBalance::endFace(SplitRound& round, int edge, int face, int& middle)
{
    if (round.done(face)) {
        return;
    }
    const Mesh& mesh = round.mesh();
    const Adjacency& adjacency = mesh.adjacency();
    const Topology& topology = mesh.topology();
    const int c = cornerIn(adjacency, edge, face);
    const int after = adjacency.next(c);
    if (middle == -1) {
        middle = middlePoint(round.points(), adjacency.vertex(c), adjacency.vertex(after));
    }
    std::vector<int> polygon;
    for (int k = topology.faceBegin(face); k < topology.faceEnd(face); ++k) {
        if (k == after) {
            polygon.push_back(middle);
        }
        polygon.push_back(topology.vertex(k));
    }
    NgonPolygonList list(polygon);
    round.queueSplit(face, list);
}

void ParityBalance::innerFace(SplitRound& round, int entry, int exit, int face, int entryMiddle, int& exitMiddle)
{
    if (round.done(face)) {
        return;
    }
    const Mesh& mesh = round.mesh();
    const Adjacency& adjacency = mesh.adjacency();
    NgonNewPoints& points = round.points();
    const int c1 = cornerIn(adjacency, entry, face);
    const int n1 = adjacency.next(c1);
    const int c2 = cornerIn(adjacency, exit, face);
    const int n2 = adjacency.next(c2);
    const int a1 = adjacency.vertex(c1);
    const int b1 = adjacency.vertex(n1);
    const int a2 = adjacency.vertex(c2);
    const int b2 = adjacency.vertex(n2);
    const int n = mesh.topology().faceSize(face);
    const int m1 = entryMiddle; // always handed on by the face before
    if (exitMiddle == -1) {
        exitMiddle = middlePoint(points, a2, b2);
    }
    const int m2 = exitMiddle;
    std::vector<std::vector<int>> polygons;
    if (c1 == n2 || n1 == c2) {
        // Entry and exit edge meet at a corner: a point inside joins both middles.
        const bool atEntryStart = c1 == n2;
        const int shared = atEntryStart ? a1 : b1;
        const int u = atEntryStart ? b1 : a1;
        const int w = atEntryStart ? a2 : b2;
        const Vec3f sum = points.position(shared) * 0.5f + (points.position(u) + points.position(w));
        const int centre = points.add({sum.x / innerDivisor, sum.y / innerDivisor, sum.z / innerDivisor});
        points.addParent(u, innerShare);
        points.addParent(w, innerShare);
        points.addParent(shared, cornerShare);
        points.endParents();
        if (atEntryStart) {
            const int y = adjacency.next(n1);
            polygons.push_back({shared, m1, centre, m2});
            polygons.push_back({centre, m1, b1, adjacency.vertex(y)});
            polygons.push_back({centre});
            cornerRun(adjacency, y, c2, polygons.back());
            polygons.back().push_back(m2);
        } else {
            const int y = adjacency.next(n2);
            polygons.push_back({shared, m2, centre, m1});
            polygons.push_back({centre, m2, b2, adjacency.vertex(y)});
            polygons.push_back({centre});
            cornerRun(adjacency, y, c1, polygons.back());
            polygons.back().push_back(m1);
        }
    } else if (n == 4) {
        polygons.push_back({m1, m2});
        cornerRun(adjacency, n2, c1, polygons.back());
        polygons.push_back({m2, m1});
        cornerRun(adjacency, n1, c2, polygons.back());
    } else {
        polygons.push_back({m1});
        cornerRun(adjacency, n1, c2, polygons.back());
        polygons.back().push_back(m2);
        cornerRun(adjacency, n2, c1, polygons.back());
    }
    round.queuePolygons(face, polygons);
}

} // namespace areko::remesher::extract
