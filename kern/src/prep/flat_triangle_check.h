#pragma once

#include "mesh/bit_array.h"
#include "mesh/mesh.h"

#include <bit>

namespace areko::remesher::prep {

// Parameters and results of one run of the flat triangle check.
struct FlatTriangleBlock {
    int mode = 3;                                          // 0 .. 3
    int round = 0;
    float angleSmall = std::bit_cast<float>(0x3c8efa35u);  // 1 degree
    float angleMain = std::bit_cast<float>(0x3eb2b8c2u);   // 20 degrees
    int requiredBegin = -1;                                // point range of required edges (mode 1)
    int requiredEnd = -1;
    // results
    int count = 0;       // number of chosen edges
    bool again = false;  // another run is worthwhile (modes 0 and 3)
};

// Flips edges at acute corners of flat triangles of a triangle mesh; edges set
// in `locked` are never flipped. In mode 1 with a point range, edges folded by
// at least 100 degrees between two points of the range are flipped first and
// without a check. Points do not move.
class FlatTriangleCheck {
public:
    static void run(mesh::Mesh& mesh, const mesh::BitArray& locked, FlatTriangleBlock& block);
};

} // namespace areko::remesher::prep
