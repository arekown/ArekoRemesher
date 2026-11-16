#include "sharp_edges.h"

#include "attribute_maps.h"
#include "link_path.h"
#include "mesh.h"


namespace areko::remesher::mesh {

namespace {

// Edge bits -> one bit per corner: every corner of a set edge.
BitArray cornersOfEdges(const Mesh& mesh, const BitArray& edges)
{
    const Adjacency& adjacency = mesh.adjacency();
    BitArray corners(static_cast<std::size_t>(mesh.cornerCount()));
    const std::size_t n = std::min(edges.size(), static_cast<std::size_t>(adjacency.edgeCount()));
    for (std::size_t e = 0; e < n; ++e) {
        if (edges.test(e)) {
            for (const int c : adjacency.edgeCorners(static_cast<int>(e))) {
                corners.set(static_cast<std::size_t>(c));
            }
        }
    }
    return corners;
}

// Corner bits -> edge bits: an edge is set when one of its corners is.
BitArray edgesOfCorners(const Mesh& mesh, const BitArray& corners)
{
    const Adjacency& adjacency = mesh.adjacency();
    BitArray bits(static_cast<std::size_t>(adjacency.edgeCount()));
    const std::size_t n = std::min(corners.size(), static_cast<std::size_t>(mesh.cornerCount()));
    for (std::size_t c = 0; c < n; ++c) {
        const int e = adjacency.cornerEdge(static_cast<int>(c));
        if (e >= 0 && corners.test(c)) {
            bits.set(static_cast<std::size_t>(e));
        }
    }
    return bits;
}

} // namespace

BoolMap* SharpEdges::edgeForm(Mesh& mesh)
{
    if (auto* edges = mesh.maps().findAs<BoolMap>(map_names::sharpEdges, ElementKind::Edge)) {
        return edges;
    }
    if (toEdgeForm(mesh)) {
        return mesh.maps().findAs<BoolMap>(map_names::sharpEdges, ElementKind::Edge);
    }
    return nullptr;
}

bool SharpEdges::toHalfEdgeForm(Mesh& mesh)
{
    const auto* edges = mesh.maps().findAs<BoolMap>(map_names::sharpEdges, ElementKind::Edge);
    if (edges == nullptr) {
        return false;
    }
    BitArray corners = cornersOfEdges(mesh, edges->bits());
    mesh.maps().remove(map_names::sharpEdges, ElementKind::Edge);
    mesh.maps().set(std::make_unique<BoolMap>(std::string(map_names::sharpEdges), ElementKind::HalfEdge, std::move(corners)));
    return true;
}

bool SharpEdges::toEdgeForm(Mesh& mesh)
{
    const auto* halfEdges = mesh.maps().findAs<BoolMap>(map_names::sharpEdges, ElementKind::HalfEdge);
    if (halfEdges == nullptr) {
        return false;
    }
    setEdgeBits(mesh, edgesOfCorners(mesh, halfEdges->bits()));
    return true;
}

void SharpEdges::setEdgeBits(Mesh& mesh, BitArray bits)
{
    mesh.maps().remove(map_names::sharpEdges, ElementKind::HalfEdge);
    mesh.maps().set(std::make_unique<BoolMap>(std::string(map_names::sharpEdges), ElementKind::Edge, std::move(bits)));
}

std::optional<BitArray> SharpEdges::halfEdgeBits(const Mesh& mesh)
{
    const AttributeMaps& maps = mesh.maps();
    if (const auto* halfEdges = maps.findAs<BoolMap>(map_names::sharpEdges, ElementKind::HalfEdge)) {
        return halfEdges->bits();
    }
    if (const auto* edges = maps.findAs<BoolMap>(map_names::sharpEdges, ElementKind::Edge)) {
        return cornersOfEdges(mesh, edges->bits());
    }
    return std::nullopt;
}

std::optional<BitArray> SharpEdges::edgeBits(const Mesh& mesh)
{
    const AttributeMaps& maps = mesh.maps();
    if (const auto* edges = maps.findAs<BoolMap>(map_names::sharpEdges, ElementKind::Edge)) {
        BitArray bits = edges->bits();
        bits.resize(static_cast<std::size_t>(mesh.adjacency().edgeCount()));
        return bits;
    }
    if (const auto* halfEdges = maps.findAs<BoolMap>(map_names::sharpEdges, ElementKind::HalfEdge)) {
        return edgesOfCorners(mesh, halfEdges->bits());
    }
    return std::nullopt;
}

BitArray SharpEdges::extendableEdgeBits(Mesh& mesh)
{
    const auto edges = static_cast<std::size_t>(mesh.adjacency().edgeCount());
    BitArray bits(edges);
    if (const BoolMap* existing = edgeForm(mesh)) {
        bits = existing->bits();
        bits.resize(edges);
    }
    return bits;
}

void SharpEdges::fromChains(Mesh& mesh, const Csr& chains, const IndexMap& points, int oldPointCount)
{
    const int pointCount = mesh.pointCount();
    std::vector<int> renumber(static_cast<std::size_t>(oldPointCount), -1);
    const std::uint8_t state = points.state & 3;
    if (state == 1 || state == 2) {
        for (int v = 0; v < oldPointCount && v < pointCount; ++v) {
            renumber[static_cast<std::size_t>(v)] = v;
        }
    } else if (points.newToOld) {
        const auto& newToOld = *points.newToOld;
        for (std::size_t i = 0; i < newToOld.size(); ++i) {
            const int o = newToOld[i];
            if (o >= 0 && o < oldPointCount && renumber[static_cast<std::size_t>(o)] < 0) {
                renumber[static_cast<std::size_t>(o)] = static_cast<int>(i);
            }
        }
    } else {
        return;
    }

    const Adjacency& adjacency = mesh.adjacency();
    const auto& positions = mesh.points();
    BitArray bits(static_cast<std::size_t>(adjacency.edgeCount()));
    auto link = [&](int a, int b) {
        const LinkPath path = shortestLink(adjacency, positions, a, b);
        if (path.first >= 0) {
            bits.set(static_cast<std::size_t>(path.first));
        }
        if (path.second >= 0) {
            bits.set(static_cast<std::size_t>(path.second));
        }
    };
    for (int r = 0; r < chains.rows(); ++r) {
        int previous = -1;
        for (const int p : chains.row(r)) {
            const int q = (p >= 0 && p < oldPointCount) ? renumber[static_cast<std::size_t>(p)] : -1;
            if (q < 0) {
                continue;
            }
            if (previous >= 0) {
                link(previous, q);
            }
            previous = q;
        }
    }
    setEdgeBits(mesh, std::move(bits));
}

} // namespace areko::remesher::mesh
