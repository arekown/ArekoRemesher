#pragma once

#include "adjacency.h"
#include "vec3.h"

#include <span>

namespace areko::remesher::mesh {

// Connection of two points over at most two edges:
//  - the first edge in pointEdges(a) whose other end is b, else
//  - over the common neighbour c with the smallest |c - a| + |c - b|, taken in
//    the edge order of b (only strictly smaller sums replace the first one).
// `first` runs from a (to b or c), `second` from c to b; -1 where absent.
struct LinkPath {
    int first = -1;
    int second = -1;
    int via = -1;

    [[nodiscard]] bool found() const noexcept { return first >= 0; }
};

[[nodiscard]] LinkPath shortestLink(const Adjacency& adjacency, std::span<const Vec3f> points, int a, int b);

} // namespace areko::remesher::mesh
