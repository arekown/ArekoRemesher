#pragma once

#include "param/rounding_step.h"

namespace areko::remesher::param {

// Rounds everything that is left: every position of the integer list that is
// not fixed yet gets the whole number nearest to its value in the latest
// solution (half away from zero), in ascending order of the positions. Grid
// periods do not count.
class FinalRounding final : public RoundingStep {
public:
    Outcome select(const SolveState& state, solver::DriverRequest& request) override;
};

} // namespace areko::remesher::param
