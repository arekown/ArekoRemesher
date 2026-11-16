#pragma once

#include "mesh/attribute_map.h"
#include "mesh/mesh.h"

#include <vector>

namespace areko::remesher::extract {

// Puts the output points of the kept border points back onto them: for every
// border loop of the triangle mesh (in loop order) and every point v of it
// with freeze mark 1 and an output point, that point takes the position of
// v, bit for bit.
class FrozenBorderReset {
public:
    static void run(mesh::Mesh& output, const mesh::Mesh& triangles, const std::vector<mesh::FreezeMark>& marks,
                    const std::vector<int>& originalToOutput);
};

} // namespace areko::remesher::extract
