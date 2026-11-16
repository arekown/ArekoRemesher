#pragma once

#include "adjacency.h"
#include "csr.h"

#include <span>

namespace areko::remesher::mesh {

// Groups of connected edges from a list. The first edge of the list that is
// not grouped yet starts a group with its end point [0] on a stack (last in,
// first out, each point at most once at a time). The point taken from the
// stack adds its still free list edges in pointEdges order to the group and
// puts their other end points on the stack.
class EdgeGroups {
public:
    // Row g lists the edges of group g in the order they joined it.
    [[nodiscard]] static Csr byConnectivity(const Adjacency& adjacency, std::span<const int> edges);
};

} // namespace areko::remesher::mesh
