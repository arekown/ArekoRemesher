#include "mesh_repair.h"

#include "mesh/boundary_loops.h"
#include "mesh/degenerate_cleanup.h"
#include "mesh/edge_groups.h"
#include "mesh/modification.h"
#include "util/introsort.h"

#include <algorithm>
#include <numeric>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

// Fan number per corner of v (in pointCorners order) and the number of fans.
// Corners belong to the same fan when an edge at v with exactly two corners
// joins them; fans are numbered by their smallest corner.
int fansOf(const Adjacency& adjacency, int v, std::vector<int>& fan)
{
    const auto corners = adjacency.pointCorners(v);
    fan.assign(corners.size(), -1);
    auto indexOf = [&](int c) {
        return static_cast<std::size_t>(std::lower_bound(corners.begin(), corners.end(), c) - corners.begin());
    };
    int count = 0;
    std::vector<int> stack;
    for (std::size_t k = 0; k < corners.size(); ++k) {
        if (fan[k] >= 0) {
            continue;
        }
        fan[k] = count;
        stack.assign(1, corners[k]);
        while (!stack.empty()) {
            const int c = stack.back();
            stack.pop_back();
            for (const int half : {c, adjacency.prev(c)}) {
                const int e = adjacency.cornerEdge(half);
                if (e < 0 || adjacency.edgeCorners(e).size() != 2) {
                    continue;
                }
                const auto pair = adjacency.edgeCorners(e);
                const int other = pair[0] == half ? pair[1] : pair[0];
                const int atV = adjacency.vertex(other) == v ? other : adjacency.next(other);
                const std::size_t i = indexOf(atV);
                if (i < corners.size() && corners[i] == atV && fan[i] < 0) {
                    fan[i] = count;
                    stack.push_back(atV);
                }
            }
        }
        ++count;
    }
    return count;
}

} // namespace

int MeshRepair::run(Mesh& mesh)
{
    int split = splitNonManifold(mesh);
    (void)orient(mesh);
    const int welded = weldLoopPoints(mesh);
    if (!mesh.adjacency().nonManifoldEdges().empty()) {
        split += splitNonManifold(mesh);
    }
    const int removed = removeUnusedPoints(mesh);
    return removed + split + welded;
}

int MeshRepair::splitNonManifold(Mesh& mesh)
{
    return static_cast<int>(splitNonManifoldPoints(mesh, nullptr).size());
}

std::vector<int> MeshRepair::splitNonManifoldPoints(Mesh& mesh, const BitArray* candidates)
{
    const Adjacency& adjacency = mesh.adjacency();
    const auto chosen = [candidates](std::size_t v) { return candidates == nullptr || candidates->testBounded(v); };
    BitArray mask = candidates != nullptr ? *candidates : BitArray(static_cast<std::size_t>(mesh.pointCount()), true);
    // At an isolated non-manifold edge the end point with fewer corners waits
    // for the second pass (the second end point on ties).
    const Csr groups = EdgeGroups::byConnectivity(adjacency, adjacency.nonManifoldEdges());
    for (int g = 0; g < groups.rows(); ++g) {
        if (groups.rowSize(g) != 1) {
            continue;
        }
        const auto ends = adjacency.edgeVertices(groups.row(g)[0]);
        const bool firstFewer = adjacency.pointCorners(ends[0]).size() < adjacency.pointCorners(ends[1]).size();
        mask.setBounded(static_cast<std::size_t>(firstFewer ? ends[0] : ends[1]), false);
    }
    std::vector<int> sources;
    splitFans(mesh, mask, sources);
    BitArray inverse(mask.size());
    for (std::size_t v = 0; v < mask.size(); ++v) {
        inverse.set(v, !mask.test(v) && chosen(v));
    }
    splitFans(mesh, inverse, sources);
    return sources;
}

void MeshRepair::splitFans(Mesh& mesh, const BitArray& mask, std::vector<int>& sources)
{
    const Adjacency& adjacency = mesh.adjacency();
    const int pointCount = mesh.pointCount();
    std::vector<int> corners = mesh.topology().cornerVertex();
    std::vector<int> added;
    std::vector<int> fan;
    const int masked = std::min(pointCount, static_cast<int>(mask.size()));
    for (int v = 0; v < masked; ++v) {
        if (!mask.test(static_cast<std::size_t>(v)) || fansOf(adjacency, v, fan) < 2) {
            continue;
        }
        const int fans = *std::max_element(fan.begin(), fan.end()) + 1;
        const int first = pointCount + static_cast<int>(added.size());
        const auto vertexCorners = adjacency.pointCorners(v);
        for (std::size_t k = 0; k < vertexCorners.size(); ++k) {
            if (fan[k] > 0) {
                corners[static_cast<std::size_t>(vertexCorners[k])] = first + fan[k] - 1;
            }
        }
        added.insert(added.end(), static_cast<std::size_t>(fans - 1), v);
    }
    if (added.empty()) {
        return;
    }
    sources.insert(sources.end(), added.begin(), added.end());
    splitPoints(mesh, std::move(corners), added);
}

void MeshRepair::splitPoints(Mesh& mesh, std::vector<int> cornerVertex, const std::vector<int>& sources)
{
    Modification change;
    Csr copySources;
    for (const int v : sources) {
        change.appendPoints.push_back(mesh.points()[static_cast<std::size_t>(v)]);
        copySources.entries.push_back(v);
        copySources.starts.push_back(static_cast<int>(copySources.entries.size()));
    }
    change.pointSources = std::move(copySources);
    change.cornerVertex = std::move(cornerVertex);
    change.faces.state = 1;
    change.corners.state = 1;
    change.halfEdges.state = 1;
    MeshModifier::apply(mesh, change);
}

bool MeshRepair::orient(Mesh& mesh)
{
    const Topology& topology = mesh.topology();
    const Adjacency& adjacency = mesh.adjacency();
    const ConnectedParts& parts = mesh.connectedParts();
    constexpr char unvisited = 2;
    std::vector<char> state(static_cast<std::size_t>(topology.faceCount()), unvisited);
    std::vector<int> queue;
    for (int t = 0; t < parts.partCount(); ++t) {
        const int start = parts.partFaces().row(t)[0];
        state[static_cast<std::size_t>(start)] = 0;
        queue.assign(1, start);
        for (std::size_t head = 0; head < queue.size(); ++head) {
            const int f = queue[head];
            const char own = state[static_cast<std::size_t>(f)];
            for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
                const int e = adjacency.cornerEdge(c);
                if (e < 0 || adjacency.edgeCorners(e).size() != 2) {
                    continue;
                }
                const auto pair = adjacency.edgeCorners(e);
                const int g = adjacency.cornerFace(pair[0] == c ? pair[1] : pair[0]);
                // Both half edges starting at the same point run the same way.
                const char expected =
                    adjacency.vertex(pair[0]) == adjacency.vertex(pair[1]) ? static_cast<char>(1 - own) : own;
                char& other = state[static_cast<std::size_t>(g)];
                if (other == unvisited) {
                    other = expected;
                    queue.push_back(g);
                } else if (other != expected) {
                    return false;
                }
            }
        }
    }
    if (std::find(state.begin(), state.end(), char{1}) == state.end()) {
        return true;
    }
    std::vector<int> corners = topology.cornerVertex();
    std::vector<int> cornerOrigin(corners.size());
    std::iota(cornerOrigin.begin(), cornerOrigin.end(), 0);
    std::vector<int> halfEdgeOrigin = cornerOrigin;
    std::vector<int> faceOrigin(state.size());
    std::iota(faceOrigin.begin(), faceOrigin.end(), 0);
    for (int f = 0; f < topology.faceCount(); ++f) {
        if (state[static_cast<std::size_t>(f)] != 1) {
            continue;
        }
        // (p0, ..., pn-1) -> (pn-1, ..., p0); the new half edge i runs back along the old half edge n-2-i.
        const int b = topology.faceBegin(f);
        const int n = topology.faceSize(f);
        for (int i = 0; i < n; ++i) {
            const auto c = static_cast<std::size_t>(b + i);
            corners[c] = topology.vertex(b + n - 1 - i);
            cornerOrigin[c] = b + n - 1 - i;
            halfEdgeOrigin[c] = b + (2 * n - 2 - i) % n;
        }
    }
    Modification change;
    change.cornerVertex = std::move(corners);
    change.points.state = 1;
    change.faces.newToOld = std::move(faceOrigin); // a renumbering also drops the UV corners, which no longer fit
    change.corners.newToOld = std::move(cornerOrigin);
    change.halfEdges.newToOld = std::move(halfEdgeOrigin);
    MeshModifier::apply(mesh, change);
    return true;
}

int MeshRepair::weldLoopPoints(Mesh& mesh)
{
    const Adjacency& adjacency = mesh.adjacency();
    const auto& points = mesh.points();
    const BoundaryLoops loops(mesh);
    std::vector<int> corners = mesh.topology().cornerVertex();
    int welded = 0;
    std::vector<Vec3f> q;
    std::vector<int> order;
    std::vector<int> representative;
    std::vector<int> fan;
    for (int r = 0; r < loops.count(); ++r) {
        const auto loop = loops.points().row(r);
        const auto n = static_cast<int>(loop.size());
        q.clear();
        for (const int v : loop) {
            q.push_back(points[static_cast<std::size_t>(v)]);
        }
        order.resize(static_cast<std::size_t>(n));
        std::iota(order.begin(), order.end(), 0);
        util::Introsort::sort(order.begin(), order.end(), [&](int a, int b) {
            return lessXyz(q[static_cast<std::size_t>(a)], q[static_cast<std::size_t>(b)]);
        });
        representative.resize(static_cast<std::size_t>(n));
        std::iota(representative.begin(), representative.end(), 0);
        bool allEqual = true;
        bool anyEqual = false;
        for (int j = 0; j + 1 < n; ++j) {
            const int a = order[static_cast<std::size_t>(j)];
            const int b = order[static_cast<std::size_t>(j) + 1];
            if (q[static_cast<std::size_t>(a)] == q[static_cast<std::size_t>(b)]) {
                representative[static_cast<std::size_t>(b)] = representative[static_cast<std::size_t>(a)];
                anyEqual = true;
            } else {
                allEqual = false;
            }
        }
        bool redirect = allEqual;
        if (!allEqual && anyEqual) {
            // Test mesh: the loop as one polygon over its representatives.
            Mesh test;
            std::vector<Vec3f> testPoints = q;
            test.swapPoints(testPoints);
            std::vector<int> start{0, n};
            std::vector<int> face = representative;
            test.editTopology().swapPolygons(start, face);
            test.topologyChanged();
            const Adjacency& testAdjacency = test.adjacency();
            redirect = testAdjacency.nonManifoldEdges().empty();
            for (int v = 0; redirect && v < n; ++v) {
                redirect = fansOf(testAdjacency, v, fan) < 2;
            }
        }
        if (!redirect) {
            continue;
        }
        for (int i = 1; i < n; ++i) {
            const int v = loop[static_cast<std::size_t>(i)];
            const int t = loop[static_cast<std::size_t>(representative[static_cast<std::size_t>(i)])];
            if (v == t) {
                continue;
            }
            ++welded;
            for (const int c : adjacency.pointCorners(v)) {
                corners[static_cast<std::size_t>(c)] = t;
            }
        }
    }
    if (welded == 0) {
        return 0;
    }
    Modification change;
    change.cornerVertex = std::move(corners);
    change.points.state = 1;
    change.faces.state = 1;
    change.corners.state = 1;
    change.halfEdges.state = 1;
    MeshModifier::apply(mesh, change);
    DegenerateCleanup::run(mesh, {.keepMaps = true});
    return welded;
}

int MeshRepair::removeUnusedPoints(Mesh& mesh, std::vector<int>* renumbering)
{
    const Adjacency& adjacency = mesh.adjacency();
    const int pointCount = mesh.pointCount();
    std::vector<int> renumber(static_cast<std::size_t>(pointCount), -1);
    std::vector<Vec3f> kept;
    for (int v = 0; v < pointCount; ++v) {
        if (!adjacency.pointCorners(v).empty()) {
            renumber[static_cast<std::size_t>(v)] = static_cast<int>(kept.size());
            kept.push_back(mesh.points()[static_cast<std::size_t>(v)]);
        }
    }
    const int removed = pointCount - static_cast<int>(kept.size());
    if (removed == 0) {
        return 0;
    }
    std::vector<int> corners = mesh.topology().cornerVertex();
    for (int& v : corners) {
        v = renumber[static_cast<std::size_t>(v)];
    }
    if (renumbering != nullptr) {
        *renumbering = renumber;
    }
    Modification change;
    change.points.count = static_cast<int>(kept.size());
    change.points.oldToNew = std::move(renumber);
    change.replacePoints = std::move(kept);
    change.cornerVertex = std::move(corners);
    change.faces.state = 1;
    change.corners.state = 1;
    change.halfEdges.state = 1;
    change.edges.state = 1;
    MeshModifier::apply(mesh, change);
    return removed;
}

} // namespace areko::remesher::prep
