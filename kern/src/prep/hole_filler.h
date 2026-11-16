#pragma once

#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Closes small border loops with one face each. A loop with 3 .. maxPoints
// points is closed when its perimeter (float sum of the edge lengths) is
// below `maxPerimeter` or its area measure 0.25 * sum |(p_i - p_i-1) x
// (p_i+1 - p_i)| is below `maxArea` (each rule only when its limit is > 0).
// The new face runs along the loop and takes the face maps of the first face
// at the first loop point; the mesh becomes a polygon mesh. Afterwards the
// points of the new faces are smoothed (list order, uniform weight, no fold
// guard).
class HoleFiller {
public:
    struct Options {
        float maxPerimeter = 0.0f;
        float maxArea = 0.0f;
        int maxPoints = 4;
        int smoothIterations = 2;
        float smoothWeight = 0.2f;
    };

    // Returns the number of closed loops.
    static int fill(mesh::Mesh& mesh, const Options& options);
};

} // namespace areko::remesher::prep
