#pragma once

#include "param/constraint_set.h"
#include "param/region_bridges.h"
#include "param/region_growth.h"
#include "param/soft_system.h"
#include "param/system_context.h"
#include "param/uv_output.h"

#include <span>

namespace areko::remesher::param {

// The rows the bridges give the additional system. A bridge compares one
// coordinate at the sources of two neighbouring regions; how far they lie
// apart along and across that coordinate decides what is asked of them:
//  - a soft row "difference = present difference", weighted,
//  - with it a hard row "difference = count of the bridge", where the count
//    is the column of the bridge and becomes a whole number of at least 1,
//  - or only a hard row "difference = 1" for two curves closer than that.
// Two curves always get weight 1; two singular points lying side by side
// across the coordinate are held together strongly, and less the further
// apart or the more askew they lie. The side on which the difference counts
// as positive follows the combed direction along the way between the
// sources. Terms of anchored variables are left out. A bridge without edge
// and, on a symmetric part, one with a source in the mirrored half gets no
// row.
class BridgeRows {
public:
    // Files the rows of all bridges, the smallest span first. `growth` holds
    // the regions of the first and of the second coordinate, `constraints`
    // the conditions of the main system.
    static void add(const SystemContext& context, const ConstraintSet& constraints, const UvOutput::Coordinates& uv,
                    std::span<const RegionGrowth, 2> growth, const RegionBridges& bridges, SoftSystem& system);
};

} // namespace areko::remesher::param
