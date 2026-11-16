#pragma once

#include "param/system.h"
#include "param/system_context.h"

#include <span>

namespace areko::remesher::param {

// Numbering of the unknowns of a system: the variables that are left after
// the conditions, in ascending order.
class UnknownNumbering {
public:
    // Writes variable -> unknown and back. A variable gets no unknown when
    // it has a condition or belongs to a copy point of the mirrored half.
    static void assign(const SystemContext& context, System& system);

    // Rewrites every term of the rows from its variable to its unknown; rows
    // sorted by variable stay sorted. Returns the number of terms whose
    // variable has no unknown (they get -1).
    [[nodiscard]] static int renumber(const System& system, std::span<TermList> rows);

    // Writes the unknowns that have to become whole numbers (every unknown
    // whose variable is not continuous), and their periods if the conditions
    // changed one.
    static void listIntegers(System& system);
};

} // namespace areko::remesher::param
