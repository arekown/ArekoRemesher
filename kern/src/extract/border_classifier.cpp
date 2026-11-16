#include "border_classifier.h"

#include "edge_dissolve.h"
#include "node_search.h"
#include "trace_node.h"

#include "mesh/boundary_loops.h"
#include "parts/part.h"

#include <algorithm>
#include <cmath>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

// Distance below which a point counts as lying on a border edge (absolute).
constexpr float onEdgeLimit = 0.0001f;

} // namespace

bool BorderClassifier::anchor(const RawMesh& raw, int point)
{
    return raw.points.onBorderEdge.testBounded(static_cast<std::size_t>(point)) &&
           raw.points.source[static_cast<std::size_t>(point)].type == node_kind::point;
}

std::vector<int> BorderClassifier::visitOrder(std::span<const int> loop, const RawMesh& raw)
{
    const auto first = std::find_if(loop.begin(), loop.end(), [&](int q) { return anchor(raw, q); });
    if (first == loop.end()) {
        return {loop.begin(), loop.end()};
    }
    std::vector<int> order(first, loop.end());
    order.insert(order.end(), loop.begin(), first + 1);
    return order;
}

void BorderClassifier::run(RawMesh& raw, const parts::Part& part)
{
    const Mesh& triangles = part.mesh;
    const Adjacency& original = triangles.adjacency();
    BitArray& onBorder = raw.points.onBorderEdge;
    const std::vector<PointSource>& source = raw.points.source;
    const BoundaryLoops loops(raw.mesh);

    std::vector<int> stretch;
    int unmarked = 0;
    for (int r = 0; r < loops.count(); ++r) {
        int last = -1;
        for (const int q : visitOrder(loops.points().row(r), raw)) {
            if (!anchor(raw, q)) {
                stretch.push_back(q);
                unmarked += onBorder.testBounded(static_cast<std::size_t>(q)) ? 0 : 1;
                continue;
            }
            if (last != -1) {
                const int a = source[static_cast<std::size_t>(last)].id;
                const int b = source[static_cast<std::size_t>(q)].id;
                const int e = original.edgeBetween(a, b);
                if (e != -1 && original.edgeCorners(e).size() == 1 && !stretch.empty() && unmarked != 0) {
                    const Vec3f pa = triangles.points()[static_cast<std::size_t>(a)];
                    const Vec3f pb = triangles.points()[static_cast<std::size_t>(b)];
                    for (const int p : stretch) {
                        const auto bit = static_cast<std::size_t>(p);
                        if (onBorder.testBounded(bit)) {
                            continue;
                        }
                        const Vec3f at = raw.mesh.points()[bit];
                        // Every point of a border loop lies in a face, so it
                        // is a node below the bit count; the extra points of
                        // a frozen border (from the bit count on) lie in no
                        // face and never come here.
                        if (std::sqrt(NodeSearch::segmentDistanceSq(at, pa, pb)) < onEdgeLimit) {
                            onBorder.setBounded(bit);
                        }
                    }
                }
            }
            last = q;
            unmarked = 0;
            stretch.clear();
        }
    }

    if (part.frozenBorderMarks() == nullptr) {
        return;
    }
    std::vector<int> edges;
    {
        const Adjacency& adjacency = raw.mesh.adjacency();
        BitArray chosen(static_cast<std::size_t>(adjacency.edgeCount()));
        for (const int q : loops.points().entries) {
            const auto bit = static_cast<std::size_t>(q);
            const PointSource from = source[bit];
            if (!onBorder.testBounded(bit) || from.type == node_kind::point) {
                continue;
            }
            if (from.type == node_kind::edge && part.lockedBorderEdges.testBounded(static_cast<std::size_t>(from.id))) {
                continue;
            }
            for (const int e : adjacency.pointEdges(q)) {
                if (adjacency.edgeCorners(e).size() > 1) {
                    chosen.set(static_cast<std::size_t>(e));
                }
            }
        }
        chosen.forEachSet([&](std::size_t e) { edges.push_back(static_cast<int>(e)); });
    }
    EdgeDissolve::run(raw.mesh, edges, true);
}

} // namespace areko::remesher::extract
