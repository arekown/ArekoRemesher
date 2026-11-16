#pragma once

#include "param/rounding_step.h"
#include "parts/part.h"
#include "run/progress.h"

#include <functional>

namespace areko::remesher::param {

// Solves the system of the parameterisation of a part and turns the solution
// into grid coordinates:
//  - the system is solved; while unknowns still have to become whole numbers
//    the rounds follow SolveStrategy, each ending with a new solution,
//  - every variable gets its value from the last solution (eliminated ones
//    through their conditions), two per point of the cut copy,
//  - on a symmetric part the mirrored half takes the coordinates of the
//    computed half with u negated,
//  - a part without singular point gets its border points with whole
//    coordinates as start points of the extraction,
//  - with symmetry switched on, singular points and start points are made
//    symmetric.
class SolutionBuilder {
public:
    struct Options {
        bool symmetry = false;      // the part may be computed by halves
        bool freezeBorders = false; // the borders of the mesh are to be frozen
        // Called with the solved share (0 .. 1) of the part as it grows.
        std::function<void(float)> report;
    };
    // The rounding steps the rounds may ask.
    struct Steps {
        RoundingStep& regions;
        RoundingStep& borders;
    };

    // Writes the result to part.solution (an earlier one is dropped). The
    // part stays without solution when it has no system, the system cannot
    // be solved or a rounding step ends the rounds without a value. The
    // weights and the normal equations of part.system are left as the last
    // round used them. A stop request throws from the progress check.
    static void compute(parts::Part& part, const Options& options, const run::Progress& progress);
    // The same with rounding steps of the caller.
    static void compute(parts::Part& part, const Options& options, const run::Progress& progress, Steps steps);
};

} // namespace areko::remesher::param
