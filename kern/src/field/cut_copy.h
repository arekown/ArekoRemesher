#pragma once

#include "mesh/adjacency.h"
#include "mesh/bit_array.h"
#include "mesh/mesh.h"
#include "mesh/topology.h"

#include <memory>
#include <span>
#include <vector>

namespace areko::remesher::field {

// A triangle mesh cut open along edges: the same triangles and corners as the
// original, but every point on a cut is doubled so that the two sides of the
// cut no longer share it. The copy has its own point and edge numbers and its
// own neighbourhood.
//  - a point becomes a candidate with its second cut or border edge (cut
//    edges in list order, then the border edges that are not cut, each with
//    its first end point first),
//  - around a candidate the triangles fall into sectors, joined across every
//    edge at the point that is not cut,
//  - the sector with the smallest corner keeps the point number, every other
//    sector gets the next new number; copy point w below the original point
//    count is therefore always a copy of original point w,
//  - without cut edges nothing is doubled, borders included.
class CutCopy {
public:
    CutCopy() = default;

    [[nodiscard]] static CutCopy cut(const mesh::Mesh& mesh, std::span<const int> cutEdges,
                                     const mesh::BitArray& cutBits);

    [[nodiscard]] bool empty() const noexcept { return adjacency_ == nullptr; }
    [[nodiscard]] int pointCount() const noexcept { return static_cast<int>(origin_.size()); }
    [[nodiscard]] const mesh::Topology& topology() const noexcept { return *topology_; }
    [[nodiscard]] const mesh::Adjacency& adjacency() const noexcept { return *adjacency_; }
    // Original point of every copy point.
    [[nodiscard]] const std::vector<int>& pointOrigins() const noexcept { return origin_; }

    // Original edges under copy edge `copyEdge` (the edges of its corners in
    // the original), each once.
    [[nodiscard]] std::vector<int> originalEdges(const mesh::Adjacency& original, int copyEdge) const;
    // Copy edges over original edge `edge` (the edges of its corners in the
    // copy): one without a cut, two along a cut.
    [[nodiscard]] std::vector<int> copyEdges(const mesh::Adjacency& original, int edge) const;
    // Copy points that stand for original point `point`: the copy points of
    // its corners in corner order, each once.
    [[nodiscard]] std::vector<int> pointCopies(const mesh::Adjacency& original, int point) const;

private:
    // The edges that the corners of edge `edge` of `from` have in `to`, each
    // once; both number the same corners.
    [[nodiscard]] static std::vector<int> edgesOfCorners(const mesh::Adjacency& from, int edge,
                                                         const mesh::Adjacency& to);

    // Both live on the heap so that the neighbourhood keeps pointing at its
    // corner list when the copy is moved.
    std::unique_ptr<mesh::Topology> topology_;
    std::unique_ptr<mesh::Adjacency> adjacency_;
    std::vector<int> origin_;
};

} // namespace areko::remesher::field
