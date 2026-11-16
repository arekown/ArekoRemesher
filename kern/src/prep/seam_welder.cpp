#include "seam_welder.h"

#include "mesh/edge_groups.h"
#include "mesh/modification.h"
#include "mesh/sharp_edges.h"
#include "util/introsort.h"

#include <cmath>
#include <numeric>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

bool finite(Vec3f p)
{
    return !std::isnan(p.x) && !std::isnan(p.y) && !std::isnan(p.z);
}

} // namespace

int SeamWelder::apply(Mesh& mesh, bool asCrease)
{
    const Adjacency& adjacency = mesh.adjacency();
    const auto& points = mesh.points();
    const std::vector<int>& borderEdges = adjacency.boundaryEdges();

    std::vector<int> candidates;
    std::vector<char> seen(static_cast<std::size_t>(mesh.pointCount()), 0);
    const Csr groups = EdgeGroups::byConnectivity(adjacency, borderEdges);
    for (int g = 0; g < groups.rows(); ++g) {
        for (const int e : groups.row(g)) {
            for (const int x : adjacency.edgeVertices(e)) {
                if (!seen[static_cast<std::size_t>(x)]) {
                    seen[static_cast<std::size_t>(x)] = 1;
                    if (finite(points[static_cast<std::size_t>(x)])) {
                        candidates.push_back(x);
                    }
                }
            }
        }
    }
    const std::size_t n = candidates.size();
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    auto at = [&](int i) { return points[static_cast<std::size_t>(candidates[static_cast<std::size_t>(i)])]; };
    util::Introsort::sort(order.begin(), order.end(), [&](int a, int b) { return lessXyz(at(a), at(b)); });
    // Representative candidate per candidate: the first one of its run of equal points.
    std::vector<int> representative = order;
    bool changed = false;
    for (std::size_t j = 0; j < n; ++j) {
        representative[static_cast<std::size_t>(order[j])] = order[j];
    }
    for (std::size_t j = 1; j < n; ++j) {
        const int a = order[j - 1];
        const int b = order[j];
        if (at(a) == at(b)) {
            representative[static_cast<std::size_t>(b)] = representative[static_cast<std::size_t>(a)];
            changed = true;
        }
    }
    if (!changed) {
        return 0;
    }

    std::vector<int> borderCorners;
    if (asCrease) {
        for (const int e : borderEdges) {
            borderCorners.push_back(adjacency.edgeCorners(e)[0]);
        }
    }
    std::vector<int> corners = mesh.topology().cornerVertex();
    int moved = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const int v = candidates[i];
        const int z = candidates[static_cast<std::size_t>(representative[i])];
        if (v == z) {
            continue;
        }
        ++moved;
        for (const int c : adjacency.pointCorners(v)) {
            corners[static_cast<std::size_t>(c)] = z;
        }
    }

    // Same faces, same corners and half edges, only their points change.
    Modification change;
    change.cornerVertex = std::move(corners);
    change.points.state = 1;
    change.faces.state = 1;
    change.corners.state = 1;
    change.halfEdges.state = 1;
    MeshModifier::apply(mesh, change);

    if (asCrease) {
        const Adjacency& welded = mesh.adjacency();
        BitArray sharp = SharpEdges::extendableEdgeBits(mesh);
        for (const int c : borderCorners) {
            const int e = welded.cornerEdge(c);
            if (e >= 0 && welded.edgeCorners(e).size() != 1) {
                sharp.set(static_cast<std::size_t>(e));
            }
        }
        SharpEdges::setEdgeBits(mesh, std::move(sharp));
    }
    return moved;
}

} // namespace areko::remesher::prep
