#pragma once

#include "param/rounding_step.h"

namespace areko::remesher::param {

// Rounds the unknowns of a part by regions. Around every curve chain and
// every singular point lies the region of the points nearest to it; between
// the sources of two neighbouring regions a whole number of grid lines has
// to fit. The step sets up a small additional system over these numbers
// (BridgeRows, SeamRows, and a row per offset across a disc cut that keeps
// it where it is), eliminates its hard rows, solves it and rounds it in
// batches; the whole numbers it finds go back to the unknowns of the main
// system. It is asked in every round after the reweighting until it
// delivers.
//  - reads: the latest solution and its coordinates on the cut copy, the
//    conditions, periods and numbering of the system, the chains, the cut
//    graph, the cross field and the sizes, and from the driver the integer
//    list, its periods and the positions already fixed,
//  - writes: the selection of the round, one entry per free position of
//    the integer list whose variable is an unknown of the additional system
//    or follows from one of its conditions; all others are left to the
//    steps behind it,
//  - returns Delivered when the additional system could be solved and
//    rounded (also with an empty selection), Failed when not - an
//    additional system without unknowns counts as not solved -, and
//    AlreadyFixed when it was solved but every position it has a value for
//    is fixed already, so that the selection stays empty.
class RegionRounding final : public RoundingStep {
public:
    // Error budget of one batch when the additional system is rounded.
    static constexpr double batchBudget = 0.5;

    Outcome select(const SolveState& state, solver::DriverRequest& request) override;
};

} // namespace areko::remesher::param
