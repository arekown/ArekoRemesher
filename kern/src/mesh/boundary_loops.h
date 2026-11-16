#pragma once

#include "csr.h"

namespace areko::remesher::mesh {

class Adjacency;
class ConnectedParts;
class Mesh;

// Closed border loops, part by part. Loop r lists its points and, in parallel,
// its border edges; edge i runs from point i to point i+1 along the direction
// of the border half edges.
//  - every part starts at its smallest border edge, later loops at the next
//    still open border edge of the part,
//  - a walk turns around the end point over inner edges to the next border
//    edge; it is dropped when it meets an edge with more than two corners or
//    a corner without edge,
//  - a loop that passes a point with more than two border edges of the part
//    is cut there into sub-loops, each starting at the earlier visit of such
//    a point, in the order they close.
class BoundaryLoops {
public:
    explicit BoundaryLoops(const Mesh& mesh);
    // The loops of a neighbourhood and its parts, whatever they belong to
    // (a mesh, a cut copy).
    BoundaryLoops(const Adjacency& adjacency, const ConnectedParts& parts);

    [[nodiscard]] const Csr& points() const noexcept { return points_; }
    [[nodiscard]] const Csr& edges() const noexcept { return edges_; }
    [[nodiscard]] int count() const noexcept { return points_.rows(); }

private:
    void split(const Adjacency& adjacency, const ConnectedParts& parts, const std::vector<char>& border, int part);

    Csr points_;
    Csr edges_;
};

} // namespace areko::remesher::mesh
