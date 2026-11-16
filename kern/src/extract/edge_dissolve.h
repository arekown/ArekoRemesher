#pragma once

#include "mesh/mesh.h"

#include <bit>
#include <cstdint>
#include <span>

namespace areko::remesher::extract {

// Dissolves edges of a polygon mesh: the faces at the given edges are joined
// into groups (connected over the given edges only, numbered by their
// smallest face), every group is replaced by one new polygon per border loop.
//  - with `removeStraightPoints`, a point of these faces whose edges outside
//    the groups are exactly two and bend by at most 30 degrees (float limit)
//    is left out as well; its other faces become groups of their own,
//  - every new polygon follows its border chain in the direction of the
//    group's first face at the first chain edge, starting at the first point
//    that stays; chains with fewer than three points give no polygon,
//  - the old faces of the groups are removed and the new polygons appended
//    in group order; the face origin of a new polygon is the smallest face of
//    its group, one entry per polygon made (a chain without polygon adds
//    none). No corner origin is handed over: the mesh change lays out its
//    own for the appended polygons. Points stay (unused ones included).
// The mesh is switched to polygon mode first. Listing an edge twice changes nothing.
class EdgeDissolve {
public:
    // 30 degrees in float.
    static constexpr float bendLimit = std::bit_cast<float>(std::uint32_t{0x3f060a92u});

    static void run(mesh::Mesh& mesh, std::span<const int> edges, bool removeStraightPoints);
};

} // namespace areko::remesher::extract
