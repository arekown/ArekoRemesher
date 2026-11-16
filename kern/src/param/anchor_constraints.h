#pragma once

#include "param/system.h"
#include "param/system_context.h"

namespace areko::remesher::param {

// Anchors of a system: in every piece of the cut copy one point whose two
// coordinates are set to zero, so that the piece cannot slide.
//  - the first copy of a singular point in the piece, singular points in
//    list order,
//  - a piece without one takes the first point of its border loops, outside
//    the mirrored half, at which curve edges of both directions meet,
//  - else the first point of its border loops outside the mirrored half,
//  - else the point of the first corner of its first triangle.
class AnchorConstraints {
public:
    // Adds two conditions without terms per piece (first coordinate, then
    // second) and closes the anchors of the system.
    static void add(const SystemContext& context, System& system);
};

} // namespace areko::remesher::param
