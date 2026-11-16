#pragma once

#include "mesh/vec3.h"

#include <array>
#include <cstdint>

namespace areko::remesher::extract {

// Where a node of the grid lies on the triangle mesh.
namespace node_kind {
inline constexpr int free = 0;  // inside a triangle
inline constexpr int point = 1; // on a point of the mesh (reference = point)
inline constexpr int edge = 2;  // on an edge of the mesh (reference = edge)
} // namespace node_kind

namespace node_flag {
// Bit d (0 .. 3): the face in direction d of this node is done.
inline constexpr std::uint8_t directions = 0x0F;
inline constexpr std::uint8_t seed = 0x10;
// On the border of the mesh: a border point or a border edge.
inline constexpr std::uint8_t border = 0x20;
// Taken over from a frozen border between two nodes; not a corner of its face.
inline constexpr std::uint8_t frozenBorder = 0x40;
// A trace from this node ran into the border.
inline constexpr std::uint8_t ranIntoBorder = 0x80;
} // namespace node_flag

// What the trace found from a node in one direction: the node it reached,
// the direction it arrived in there, the triangle and grid coordinates at
// the end, and whether the trace came back to its own grid line.
struct NodeLink {
    int node = -1;
    int direction = 0;
    int triangle = -1;
    std::array<float, 2> uv{};
    bool degenerate = false;
};

// A point of the quad grid on the triangle mesh: a place with whole (u, v)
// where the traced grid lines meet. Directions 0 .. 3: u grows, v grows, u
// falls, v falls (in the frame of the node's triangle).
struct TraceNode {
    mesh::Vec3f position;
    // Triangle of the node; -1 for seed nodes, which belong to a point.
    int triangle = -1;
    int kind = node_kind::free;
    int reference = -1;
    std::array<float, 2> uv{};
    // Number of the output point (-1 until the node is part of a face) and
    // the number it was merged into (-1 if none).
    int id = -1;
    int mergeTarget = -1;
    std::uint8_t flags = 0;
    std::array<NodeLink, 4> links;
    // The next older node in the list of the same triangle (-1 at the end).
    int older = -1;

    [[nodiscard]] bool seed() const noexcept { return (flags & node_flag::seed) != 0; }
    [[nodiscard]] bool done(int direction) const noexcept { return ((flags >> direction) & 1u) != 0; }
};

} // namespace areko::remesher::extract
