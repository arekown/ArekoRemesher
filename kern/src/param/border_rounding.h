#pragma once

#include "param/rounding_step.h"

namespace areko::remesher::param {

// Rounds the coordinates along frozen borders. It is asked once, in the
// round after the regions delivered, and only for a part with frozen borders.
//  - Border loop by border loop of the part mesh, first for u, then for v,
//    it lists that coordinate of every copy of every kept border point in
//    walking order, with its value in the latest solution (an eliminated
//    coordinate as the sum its condition names) and the nearest whole number
//    (half away from zero).
//  - Two neighbours in the list (the last and the first are neighbours too)
//    that land on the same whole number and are joined by an edge of the cut
//    copy are in conflict; one of them then keeps its value as it is: of two
//    coordinates of a frozen border that are both unknowns, the one farther
//    from its whole number (the first one when both are equally far); if
//    only the first is such an unknown, the first; otherwise the second,
//    provided it is a coordinate of a frozen border.
//  - Every entry that is an unknown with a free position in the integer list
//    is selected: with its value as it is after a conflict, else with its
//    whole number. An unknown outside the integer list is left out, and a
//    position is selected once: a point in two border loops counts in the
//    first of them.
//  - reads: the latest solution, the kinds, conditions and numbering of the
//    system, the border loops and freeze marks of the part mesh, the copies
//    of its points and the edges of the cut copy, and from the driver the
//    integer list and the positions already fixed,
//  - writes: the selection of the round,
//  - returns Delivered when it selected something, Failed when not. The run
//    goes on with the final rounding in the same round when the selection is
//    empty.
class BorderRounding final : public RoundingStep {
public:
    Outcome select(const SolveState& state, solver::DriverRequest& request) override;
};

} // namespace areko::remesher::param
