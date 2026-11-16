#pragma once

#include "part.h"

namespace areko::remesher::parts {

// Brings the edge loop and crease curves of one part into its mesh:
//  - cut: every such curve is laid onto the part mesh once, all pieces are
//    tidied and cut in; the active curves take the cut chain (mesh points
//    and their positions); a real change makes the mirror pairs old,
//  - densify: points of the cut chains, of the border and (with one mirror
//    axis and a symmetric part) of the mirror plane are feature points, the
//    chain, border and plane edges feature edges; an edge between two
//    feature points that is no feature edge is split at its middle (with its
//    mirror edge), a triangle with three feature edges at its centroid; the
//    new points then move away from slivers, the old ones stay.
class PartCurveCut {
public:
    // `mirrorAxis`: the symmetry names exactly one mirror axis (no detection, no kept half).
    static void run(Part& part, bool mirrorAxis);

private:
    static void cut(Part& part);
    static void densify(Part& part, bool mirrorAxis);
};

} // namespace areko::remesher::parts
