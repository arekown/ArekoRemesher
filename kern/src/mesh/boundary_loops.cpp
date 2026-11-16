#include "boundary_loops.h"

#include "mesh.h"

namespace areko::remesher::mesh {

namespace {

// Drops the entries after the last closed row.
void dropOpen(Csr& csr)
{
    csr.entries.resize(static_cast<std::size_t>(csr.starts.back()));
}

void close(Csr& csr)
{
    csr.starts.push_back(static_cast<int>(csr.entries.size()));
}

} // namespace

BoundaryLoops::BoundaryLoops(const Mesh& mesh) : BoundaryLoops(mesh.adjacency(), mesh.connectedParts())
{
}

BoundaryLoops::BoundaryLoops(const Adjacency& adjacency, const ConnectedParts& parts)
{
    const int edgeCount = adjacency.edgeCount();
    auto partOfEdge = [&](int e) { return parts.partOfFace(adjacency.cornerFace(adjacency.edgeCorners(e)[0])); };

    std::vector<char> border(static_cast<std::size_t>(edgeCount), 0);
    std::vector<int> remaining(static_cast<std::size_t>(parts.partCount()), 0);
    std::vector<int> first(static_cast<std::size_t>(parts.partCount()), -1);
    for (int e = 0; e < edgeCount; ++e) {
        if (adjacency.edgeCorners(e).size() != 1) {
            continue;
        }
        border[static_cast<std::size_t>(e)] = 1;
        const int t = partOfEdge(e);
        ++remaining[static_cast<std::size_t>(t)];
        if (first[static_cast<std::size_t>(t)] < 0) {
            first[static_cast<std::size_t>(t)] = e;
        }
    }
    // A walk never needs more steps than there are corners.
    const std::size_t stepLimit = adjacency.cornerEdges().size() + 1;

    for (int t = 0; t < parts.partCount(); ++t) {
        std::vector<char> open = border;
        int searchFrom = 0;
        int e = first[static_cast<std::size_t>(t)];
        int& left = remaining[static_cast<std::size_t>(t)];
        auto take = [&](int edge, int corner) {
            points_.entries.push_back(adjacency.vertex(corner));
            edges_.entries.push_back(edge);
            if (open[static_cast<std::size_t>(edge)]) {
                open[static_cast<std::size_t>(edge)] = 0;
                --left;
            }
        };
        while (e >= 0) {
            const int startEdge = e;
            int corner = adjacency.edgeCorners(e)[0];
            take(e, corner);
            std::size_t steps = 0;
            while (true) {
                int c1 = adjacency.next(corner);
                int e1 = adjacency.cornerEdge(c1);
                bool lost = e1 < 0;
                while (!lost && adjacency.edgeCorners(e1).size() == 2 && ++steps < stepLimit) {
                    c1 = adjacency.next(adjacency.otherCorner(e1, c1));
                    e1 = adjacency.cornerEdge(c1);
                    if (e1 < 0) {
                        lost = true;
                        break;
                    }
                }
                if (lost || steps >= stepLimit) {
                    dropOpen(points_);
                    dropOpen(edges_);
                    break;
                }
                if (e1 == startEdge) {
                    close(points_);
                    close(edges_);
                    split(adjacency, parts, border, t);
                    break;
                }
                if (adjacency.edgeCorners(e1).size() == 1 &&
                    points_.entries.size() - static_cast<std::size_t>(points_.starts.back()) <= border.size()) {
                    take(e1, c1);
                    corner = c1;
                    continue;
                }
                dropOpen(points_);
                dropOpen(edges_);
                if (open[static_cast<std::size_t>(e1)]) {
                    open[static_cast<std::size_t>(e1)] = 0;
                    --left;
                }
                break;
            }
            e = -1;
            if (left == 0) {
                break;
            }
            for (int k = searchFrom; k < edgeCount; ++k) {
                if (open[static_cast<std::size_t>(k)] && partOfEdge(k) == t) {
                    e = k;
                    searchFrom = k + 1;
                    break;
                }
            }
        }
    }
}

void BoundaryLoops::split(const Adjacency& adjacency, const ConnectedParts& parts, const std::vector<char>& border,
                          int part)
{
    const int row = points_.rows() - 1;
    const std::vector<int> loopPoints(points_.row(row).begin(), points_.row(row).end());
    const std::vector<int> loopEdges(edges_.row(row).begin(), edges_.row(row).end());
    const int m = static_cast<int>(loopPoints.size());

    std::vector<int> clamps;  // positions of points with more than two border edges of the part
    std::vector<int> clampId; // per clamp position: number of its point among the distinct clamp points
    std::vector<int> distinct;
    for (int i = 0; i < m; ++i) {
        const int v = loopPoints[static_cast<std::size_t>(i)];
        int z = 0;
        for (const int e : adjacency.pointEdges(v)) {
            if (border[static_cast<std::size_t>(e)] &&
                parts.partOfFace(adjacency.cornerFace(adjacency.edgeCorners(e)[0])) == part) {
                ++z;
            }
        }
        if (z > 2) {
            clamps.push_back(i);
            int id = -1;
            for (std::size_t k = 0; k < distinct.size(); ++k) {
                if (distinct[k] == v) {
                    id = static_cast<int>(k);
                }
            }
            if (id < 0) {
                id = static_cast<int>(distinct.size());
                distinct.push_back(v);
            }
            clampId.push_back(id);
        }
    }
    if (clamps.empty()) {
        return;
    }
    points_.starts.pop_back();
    edges_.starts.pop_back();
    dropOpen(points_);
    dropOpen(edges_);

    const int p = static_cast<int>(clamps.size());
    auto append = [&](int from, int to) {
        for (int i = from; i < to; ++i) {
            points_.entries.push_back(loopPoints[static_cast<std::size_t>(i)]);
            edges_.entries.push_back(loopEdges[static_cast<std::size_t>(i)]);
        }
    };
    std::vector<int> lastVisit(distinct.size(), -1);
    std::vector<int> onward(static_cast<std::size_t>(p), -1);
    std::vector<char> visited(static_cast<std::size_t>(p), 0);
    for (int s = 1; s <= p + 1; ++s) {
        const int k = (s - 1) % p;
        const int x = clampId[static_cast<std::size_t>(k)];
        const int j = lastVisit[static_cast<std::size_t>(x)];
        if (j != -1 && !visited[static_cast<std::size_t>(j)]) {
            int q = j;
            while (true) {
                visited[static_cast<std::size_t>(q)] = 1;
                if (q == p - 1) {
                    append(clamps[static_cast<std::size_t>(q)], m);
                    append(0, clamps[0]);
                    break;
                }
                append(clamps[static_cast<std::size_t>(q)], clamps[static_cast<std::size_t>(q) + 1]);
                const int next = onward[static_cast<std::size_t>(q) + 1];
                if (next == -1 || next == k) {
                    break;
                }
                q = next;
            }
            onward[static_cast<std::size_t>(j)] = k;
            close(points_);
            close(edges_);
        }
        lastVisit[static_cast<std::size_t>(x)] = k;
    }
}

} // namespace areko::remesher::mesh
