#include "cut_copy.h"

#include "util/disjoint_sets.h"

#include <algorithm>
#include <numeric>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

CutCopy CutCopy::cut(const Mesh& mesh, std::span<const int> cutEdges, const BitArray& cutBits)
{
    const Adjacency& adjacency = mesh.adjacency();
    const int pointCount = mesh.pointCount();
    std::vector<int> cornerPoint = mesh.topology().cornerVertex();
    CutCopy copy;
    copy.origin_.resize(static_cast<std::size_t>(pointCount));
    std::iota(copy.origin_.begin(), copy.origin_.end(), 0);

    if (!cutEdges.empty()) {
        // Points that meet a second cut or border edge, in the order they do.
        std::vector<int> met(static_cast<std::size_t>(pointCount), 0);
        std::vector<int> candidates;
        auto count = [&](int e) {
            for (const int v : adjacency.edgeVertices(e)) {
                int& n = met[static_cast<std::size_t>(v)];
                if (n == 1) {
                    candidates.push_back(v);
                }
                ++n;
            }
        };
        for (const int e : cutEdges) {
            count(e);
        }
        for (const int e : adjacency.boundaryEdges()) {
            if (!cutBits.testBounded(static_cast<std::size_t>(e))) {
                count(e);
            }
        }

        for (const int v : candidates) {
            const auto corners = adjacency.pointCorners(v);
            // Number of the triangle of an edge corner among the corners at v.
            auto local = [&](int edgeCorner) {
                const int atPoint = adjacency.vertex(edgeCorner) == v ? edgeCorner : adjacency.next(edgeCorner);
                return static_cast<int>(std::lower_bound(corners.begin(), corners.end(), atPoint) - corners.begin());
            };
            // The triangles around v as sets, joined across the edges that are not cut.
            util::DisjointSets sectors(corners.size());
            for (const int e : adjacency.pointEdges(v)) {
                const auto edgeCorners = adjacency.edgeCorners(e);
                if (cutBits.testBounded(static_cast<std::size_t>(e)) || edgeCorners.size() < 2) {
                    continue;
                }
                const int first = local(edgeCorners[0]);
                for (std::size_t k = 1; k < edgeCorners.size(); ++k) {
                    sectors.join(first, local(edgeCorners[k]));
                }
            }
            // The sector of the smallest corner keeps the point; the others
            // get new points in the order of their smallest corner.
            std::vector<int> pointOfSector(corners.size(), -1);
            for (std::size_t j = 0; j < corners.size(); ++j) {
                const std::size_t sector = static_cast<std::size_t>(sectors.smallest(static_cast<int>(j)));
                if (sector == 0) {
                    continue;
                }
                if (pointOfSector[sector] == -1) {
                    pointOfSector[sector] = static_cast<int>(copy.origin_.size());
                    copy.origin_.push_back(v);
                }
                cornerPoint[static_cast<std::size_t>(corners[j])] = pointOfSector[sector];
            }
        }
    }

    copy.topology_ = std::make_unique<Topology>();
    copy.topology_->swapCorners(cornerPoint, mesh.topology().mode());
    copy.adjacency_ = std::make_unique<Adjacency>(*copy.topology_, copy.pointCount());
    return copy;
}

std::vector<int> CutCopy::edgesOfCorners(const Adjacency& from, int edge, const Adjacency& to)
{
    std::vector<int> out;
    for (const int c : from.edgeCorners(edge)) {
        const int e = to.cornerEdge(c);
        if (std::find(out.begin(), out.end(), e) == out.end()) {
            out.push_back(e);
        }
    }
    return out;
}

std::vector<int> CutCopy::originalEdges(const Adjacency& original, int copyEdge) const
{
    return edgesOfCorners(*adjacency_, copyEdge, original);
}

std::vector<int> CutCopy::copyEdges(const Adjacency& original, int edge) const
{
    return edgesOfCorners(original, edge, *adjacency_);
}

std::vector<int> CutCopy::pointCopies(const Adjacency& original, int point) const
{
    std::vector<int> out;
    for (const int c : original.pointCorners(point)) {
        const int w = topology_->vertex(c);
        if (std::find(out.begin(), out.end(), w) == out.end()) {
            out.push_back(w);
        }
    }
    return out;
}

} // namespace areko::remesher::field
