#pragma once

#include "parts/part.h"
#include "run/progress.h"

namespace areko::remesher::field {

// Cut graph of a part: cuts the part open into a disc and combs the cross
// field over it.
//  - a tree over the mesh edges joins the singular points of the cross field,
//    the border loops and the mirror plane; curve edges are cheap, so the
//    tree prefers to run along them,
//  - the mesh is cut open along the tree; further cuts are added until the
//    copy is a disc,
//  - with a mirrored half the cuts of the computed half are mirrored,
//  - the cross field is combed over the cut mesh, the quarter turns left
//    across every edge are noted,
//  - curve edges are assigned to the combed direction they follow; where two
//    curves at a point would take the same direction, a curve check cut and
//    a curve check source are added and the cut graph is built once more.
class CutGraphBuilder {
public:
    struct Options {
        bool symmetry = false;      // the part may be computed by halves
        bool followBorders = false; // border edges count as curve edges
    };

    // Writes the result to part.cutGraph (an earlier one is dropped) and the
    // mirror evaluation to part.mirror.classes. The part has a cross field.
    static void compute(parts::Part& part, const Options& options, const run::Progress& progress);
};

} // namespace areko::remesher::field
