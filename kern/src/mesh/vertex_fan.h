#pragma once

#include "adjacency.h"

#include <vector>

namespace areko::remesher::mesh {

// Edges around point v in fan order, starting with `startEdge` (triangle meshes).
//  - start: the corner of `startEdge` whose point is v (first corner first);
//    a border edge whose single corner does not start at v gives an empty fan,
//  - step: from corner c over the edge of prev(c) to the other corner of that
//    edge; when that corner is not at v (flipped orientation) its next corner
//    is taken and the following step uses the edge of the corner itself,
//  - the fan ends at a missing edge, at an edge without exactly two corners
//    (included), or when it reaches the start corner again (the start edge is
//    repeated as last entry),
//  - more than 2000 corners give an empty list.
[[nodiscard]] std::vector<int> vertexFanEdges(const Adjacency& adjacency, int v, int startEdge);

// The same fan with its corners: corner i lies between edges i and i + 1
// (the start corner first); with more than 2000 corners the corner list
// is followed by all corners of the point and the edge list is empty.
// Without `limited` the fan is walked to its end however many corners it has.
struct VertexFan {
    std::vector<int> corners;
    std::vector<int> edges;
};
[[nodiscard]] VertexFan vertexFan(const Adjacency& adjacency, int v, int startEdge, bool limited = true);
// Fan around the point of `corner`, starting with that corner and its edge.
// Without `forward` the fan runs the other way: the first step crosses the
// edge of `corner` itself instead of the edge of its previous corner.
[[nodiscard]] VertexFan vertexFanFrom(const Adjacency& adjacency, int corner, bool forward = true);
// Fan around v in the same order from a fixed start: on a point with as many
// corners as edges the fan from its first corner when it takes in every
// corner; otherwise the fan from the first corner is walked backwards (over
// the edge of each corner itself) and the fan is taken from the last corner
// reached, so that an open fan starts at its border.
[[nodiscard]] VertexFan orderedVertexFan(const Adjacency& adjacency, int v);

} // namespace areko::remesher::mesh
