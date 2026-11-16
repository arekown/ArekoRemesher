#pragma once

#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Removes faces that repeat another face. Points are visited in ascending
// order, the corners of a point in its corner order; for two corners of the
// point (earlier i, later j, faces not removed yet) of faces with the same
// corner count, the later face is removed when its points, read forwards or
// backwards from the common point, match the earlier face. A face whose point
// appears twice is compared with itself, so a face like (a, a, b) goes too.
// Points stay, also when they lose their last face.
class DuplicateFaces {
public:
    // Returns the number of removed faces.
    static int remove(mesh::Mesh& mesh);
};

} // namespace areko::remesher::prep
