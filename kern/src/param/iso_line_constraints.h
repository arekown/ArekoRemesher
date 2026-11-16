#pragma once

#include "param/system.h"
#include "param/system_context.h"

#include <vector>

namespace areko::remesher::param {

// Iso lines along the curve chains: all points of a chain get the same value
// of the coordinate that stays constant along it.
//  - the reference of a chain is the first of these variables that is free;
//    without a free one the first variable of the chain stands in,
//  - a free variable of the chain becomes "= reference": a new condition of
//    the system; in the end the new conditions are put into each other,
//  - for a variable that already has a condition a row "its expression minus
//    that of the reference = 0" is left for the elimination, and the
//    reference becomes a preferred variable,
//  - along a crease or border curve the variables also have to become whole
//    numbers.
class IsoLineConstraints {
public:
    struct Result {
        // Rows for the elimination and the variables that get its first columns.
        std::vector<TermList> rows;
        std::vector<int> preferred;
    };
    // Appends the new conditions to those of the system and writes the
    // reference per curve chain and the flags.
    [[nodiscard]] static Result build(const SystemContext& context, System& system);
};

} // namespace areko::remesher::param
