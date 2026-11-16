#pragma once

#include "param/soft_system.h"
#include "param/system_context.h"
#include "param/uv_output.h"

namespace areko::remesher::param {

// Hard rows of the additional system across cut chains: where the two sides
// of a chain without a turn lie almost on the same grid line over a long
// stretch, they are put on the same line.
//  - only chains whose quarter turns are a multiple of four, that have
//    sides and do not run along a disc cut,
//  - the offset between the first points of the two sides is split into its
//    larger and its smaller part; the row is for the coordinate of the
//    smaller part (the second one on ties),
//  - the smaller part has to stay below one grid step with a long offset or
//    a long chain, or below two with a very long offset and chain,
//  - the row sets that coordinate equal at the first points of both sides.
class SeamRows {
public:
    static void add(const SystemContext& context, const UvOutput::Coordinates& uv, SoftSystem& system);
};

} // namespace areko::remesher::param
