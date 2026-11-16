#pragma once

#include "parts/part.h"
#include "run/progress.h"

namespace areko::remesher::field {

// Cross field of a part: per triangle a cross of four directions that follows
// the curvature where it is trusted and turns as little as possible elsewhere.
//  - triangles with a trusted curvature direction are anchors,
//  - every triangle gets a frame; the cross is an angle in it, and across
//    every edge the crosses may differ by a whole number of quarter turns,
//  - the quarter turns along a spanning tree from the anchors are settled at
//    once, the others are rounded step by step while the angles are solved
//    again and again,
//  - points around which the quarter turns do not add up are the singular
//    points: there the quads will not meet four at a time,
//  - with a mirrored half only one half is solved; the other takes the
//    mirrored angles.
class CrossFieldSolver {
public:
    struct Options {
        bool symmetry = false; // the part may be computed by halves
    };

    // Writes the result to part.crossField; it stays empty for a part without
    // curvature or when the system cannot be solved. Marks the anchors'
    // mirror triangles in the curvature result.
    static void compute(parts::Part& part, const Options& options, const run::Progress& progress);
};

} // namespace areko::remesher::field
