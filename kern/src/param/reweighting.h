#pragma once

#include "param/system.h"
#include "param/system_context.h"
#include "run/progress.h"

#include <span>

namespace areko::remesher::param {

// Raises the weights of the energy rows of triangles that a first solution
// maps badly, so that the next solution corrects them:
//  - a triangle that is turned over in the grid gets 5 on all four rows,
//  - where a grid step is longer than 1.8 times the size of the triangle,
//    the row along that direction gets up to 3 (reached at 4 times the size),
//  - where the two grid directions are more than 15 degrees away from a
//    right angle, the two rows across get 1 to 4 (at 90 degrees),
//  - no factor goes above 10.
// The factors are then averaged over neighbouring triangles (a crease or
// border curve passes on only the rows along it; across an odd number of
// quarter turns the rows of the two directions change places), the weights
// multiplied by them and the normal equations built anew from the kept
// transposed rows.
class Reweighting {
public:
    // Rounds of averaging. An even number: the factor of a row that belongs
    // to no triangle (a row that could not be eliminated) is 0 after every
    // odd round and would wipe out its weight.
    static constexpr int smoothingPasses = 10;
    static_assert(smoothingPasses % 2 == 0);

    // `x` is the solution the triangles are judged by. Changes the weights
    // and the normal equations of the system.
    static void apply(const SystemContext& context, System& system, std::span<const double> x,
                      const run::Progress& progress);
};

} // namespace areko::remesher::param
