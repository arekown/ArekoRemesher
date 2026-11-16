#pragma once

#include "curvature.h"

#include "parts/part.h"
#include "run/progress.h"

namespace areko::remesher::field {

// Curvature of a part: per triangle the direction the quads should follow
// and how far that direction can be trusted.
//  - the curvature is measured on spheres of several radii (levels); per
//    level every point sums the edge tensors inside its sphere, every
//    triangle takes the mean of its corners and its eigen decomposition,
//  - per triangle the most stable level wins (directions that agree with the
//    neighbouring levels), with fallbacks when none is stable,
//  - guide curves, borders (when followed) and the mirror plane fix the
//    direction of their triangles,
//  - with a mirrored half only one half is computed; the other copies the
//    values (not the direction) of its mirror triangles.
class CurvatureField {
public:
    struct Options {
        bool symmetry = false;      // the mirror plane fixes directions
        bool halfComputed = false;  // one half is computed, the other mirrored
        bool followBorders = false; // border edges fix directions
        bool severalParts = false;
    };

    // Writes the result to part.curvature; it stays empty when the largest
    // radius is 0. Gives curves without a width their width and lays the
    // bands onto the part mesh.
    static void compute(parts::Part& part, const Options& options, const run::Progress& progress);
};

} // namespace areko::remesher::field
