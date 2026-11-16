#pragma once

#include "param/system.h"
#include "param/system_context.h"

#include <cstdint>
#include <vector>

namespace areko::remesher::param {

// Kinds and grid periods of the variables of a system: which of them have to
// become whole numbers, and why.
class VariableFlags {
public:
    // Gives a variable that is still continuous the kind `flag` and the
    // period 1; a variable that is marked already keeps what it has.
    static void markIfContinuous(System& system, int variable, std::uint8_t flag);

    // Sizes the system for the cut copy and the chains of the part and marks,
    // each with period 1:
    //  - both coordinates of every copy of a node that is a singular point or
    //    has a border edge (whatever they were marked before),
    //  - the offsets of the chains along the disc cuts,
    //  - the still continuous coordinates of the copies of the sources of
    //    the curve check,
    //  - with frozen borders the still continuous coordinates of the copies
    //    of the border points outside the mirrored half, except points that
    //    lie between kept ones.
    static void assign(const SystemContext& context, System& system);

    // Prepares the elimination at the chain ends and returns the variables
    // that are to get its first columns:
    //  - the continuous coordinates of all copies of the nodes,
    //  - where a chain without offset ends in a disc cut, the copy points at
    //    the two triangles of that cut edge whose coordinates are both marked
    //    become kept free (they are eliminated last); with such a point on
    //    either side of the first edge the last edge is left as it is,
    //  - every marked coordinate of the second and further copies of a node
    //    becomes a further copy, unless the copy is kept free.
    [[nodiscard]] static std::vector<int> markChainEnds(const SystemContext& context, System& system);
};

} // namespace areko::remesher::param
