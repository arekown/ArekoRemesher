#pragma once

#include "csr.h"
#include "topology.h"

#include <array>
#include <optional>
#include <span>
#include <vector>

namespace areko::remesher::mesh {

// Neighbourhood tables of a mesh. The numbering rules decide the order of
// edges and of every list below, and with it the result of later stages:
//  - corners per point ascending,
//  - edges numbered by their smaller end point ascending; for the same smaller
//    end point the first appearance in that point's corner list wins, the
//    following half edge of a corner before the preceding one,
//  - degenerate half edges (both ends equal) get no edge (-1),
//  - corners per edge and edges per point ascending.
// Boundary lists and boundary points are built on first use. The tables refer
// to the topology they were built from; its owner keeps it alive and unchanged.
class Adjacency {
public:
    Adjacency(const Topology& topology, int pointCount);
    Adjacency(const Adjacency&) = delete;
    Adjacency& operator=(const Adjacency&) = delete;

    // The owner moved the same topology to `topology`.
    void attach(const Topology& topology) noexcept { topology_ = &topology; }

    [[nodiscard]] int pointCount() const noexcept { return pointCount_; }
    [[nodiscard]] int edgeCount() const noexcept { return edgeCount_; }

    [[nodiscard]] int next(int c) const noexcept;
    [[nodiscard]] int prev(int c) const noexcept;
    [[nodiscard]] int vertex(int c) const noexcept { return topology_->vertex(c); }
    [[nodiscard]] int cornerFace(int c) const noexcept { return cornerFace_[static_cast<std::size_t>(c)]; }
    [[nodiscard]] int cornerEdge(int c) const noexcept { return cornerEdge_[static_cast<std::size_t>(c)]; }
    [[nodiscard]] std::array<int, 2> edgeVertices(int e) const noexcept
    {
        return {edgeVertices_[2 * static_cast<std::size_t>(e)], edgeVertices_[2 * static_cast<std::size_t>(e) + 1]};
    }
    [[nodiscard]] int otherEnd(int e, int v) const noexcept;
    [[nodiscard]] std::span<const int> edgeCorners(int e) const { return edgeCorners_.row(e); }
    [[nodiscard]] std::span<const int> pointCorners(int v) const { return pointCorners_.row(v); }
    [[nodiscard]] std::span<const int> pointEdges(int v) const { return pointEdges_.row(v); }
    // Other end point for every entry of pointEdges(v).
    [[nodiscard]] std::span<const int> pointNeighbours(int v) const;

    [[nodiscard]] const std::vector<int>& cornerEdges() const noexcept { return cornerEdge_; }
    [[nodiscard]] const std::vector<int>& edgeVertexList() const noexcept { return edgeVertices_; }

    // First edge in pointEdges(a) whose other end is b, or -1.
    [[nodiscard]] int edgeBetween(int a, int b) const;
    // Edge of the first corner of f (ascending) whose edge also has a corner in g, or -1.
    [[nodiscard]] int commonEdge(int f, int g) const;
    // Face across e from f: e has exactly two corners and f owns one of them
    // (the face of the other one; f when f owns both), else -1.
    [[nodiscard]] int otherFace(int e, int f) const;
    // Corner across e from corner c: e has exactly two corners and c is one
    // of them (the other one), else -1.
    [[nodiscard]] int otherCorner(int e, int c) const;
    // The faces of the two corners of an inner edge (one with exactly two
    // corners) in corner order, or none for any other edge.
    [[nodiscard]] std::optional<std::array<int, 2>> innerEdgeFaces(int e) const;

    // Edges with exactly one corner, ascending.
    [[nodiscard]] const std::vector<int>& boundaryEdges() const;
    // Edges with more than two corners, ascending.
    [[nodiscard]] const std::vector<int>& nonManifoldEdges() const;
    // End points of the boundary edges in list order, each once at its first appearance.
    [[nodiscard]] const std::vector<int>& boundaryPoints() const;
    // Points without any corner, ascending.
    [[nodiscard]] const std::vector<int>& isolatedPoints() const;
    [[nodiscard]] bool isBoundaryEdge(int e) const { return edgeCorners_.rowSize(e) == 1; }
    // First boundary edge in pointEdges(v) other than `except`, or -1.
    [[nodiscard]] int firstBoundaryEdgeAt(int v, int except = -1) const;
    // First corner of v (pointCorners order) whose outgoing half edge (to the
    // next corner) resp. incoming half edge (from the previous corner) lies on
    // a boundary edge, or -1.
    [[nodiscard]] int boundaryCornerAt(int v, bool outgoing) const;

    struct Counts {
        int quads = 0;
        int triangles = 0;
        int polygons = 0;
        int irregularInnerPoints = 0; // inner points with an edge count other than 4
    };
    [[nodiscard]] const Counts& counts() const;

private:
    void buildCornerFaces();
    void buildCornerEdges();
    void buildEdges();
    void buildPointEdges();
    void buildBoundaryLists() const;

    const Topology* topology_;
    int pointCount_;
    std::vector<int> cornerFace_;
    std::vector<int> cornerEdge_;
    std::vector<int> edgeVertices_;
    int edgeCount_ = 0;
    Csr edgeCorners_;
    Csr pointCorners_;
    Csr pointEdges_;
    std::vector<int> pointNeighbours_;

    mutable std::optional<std::vector<int>> boundaryEdges_;
    mutable std::vector<int> nonManifoldEdges_;
    mutable std::optional<std::vector<int>> boundaryPoints_;
    mutable std::optional<std::vector<int>> isolatedPoints_;
    mutable std::optional<Counts> counts_;
};

} // namespace areko::remesher::mesh
