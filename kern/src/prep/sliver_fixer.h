#pragma once

#include "mesh/bit_array.h"
#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Moves the widest corner of tiny triangles with a corner of at most 1.5 or
// at least 175 degrees by a local Laplace step (smoother with fold guard), in
// up to 4 x 4 rounds with the step 0.3, 0.15, 0.1, 0.1. A round is repeated
// while a moved triangle's normal turned by more than acos 0.9; moves are
// never taken back. Points in `lockedPoints`, points at edges without exactly
// two corners and points whose triangles do not face the same way stay put.
// Triangle normals and doubled areas follow the moves only at the start of a
// round; degenerate triangles take replacement normals at the start of every
// outer round.
class SliverFixer {
public:
    static void run(mesh::Mesh& mesh, const mesh::BitArray& lockedPoints);
};

} // namespace areko::remesher::prep
