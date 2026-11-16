#pragma once

#include "param/constraint.h"

#include <span>
#include <vector>

namespace areko::remesher::param {

// Grid periods of the variables: an integer variable with period m is rounded
// to multiples of m.
class ConstraintPeriods {
public:
    // Makes the periods fit the conditions from `first` on (-1 = none), in one
    // pass; conditions whose pivot has period 0 are left out. For a term whose
    // share period(variable) * coefficient / period(pivot) is not whole, the
    // smallest factor from 2 to 49 that makes it whole is looked for; the
    // period of the variable becomes the least common multiple of itself and
    // that factor. Returns how many periods of 1 were changed.
    [[nodiscard]] static int adjust(std::vector<int>& periods, std::span<const Constraint> constraints, int first);
};

} // namespace areko::remesher::param
