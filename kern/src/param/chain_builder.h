#pragma once

#include "parts/part.h"
#include "run/progress.h"

namespace areko::remesher::param {

// Chains of a part for its parameterisation, from its cut graph.
//  - the cut edges are chained between nodes: ends, branch points, the
//    singular points of the cross field and the border points on a cut;
//    after every chain the next ones start at the node it ended at,
//  - a part without singular point takes one node in their place (with a
//    mirrored half the first node that touches the computed half),
//  - every chain gets the points of the cut copy on its two sides and the
//    quarter turns between the sides,
//  - marked are: the nodes that are singular points; chains of one edge at a
//    singular point with a whole turn whose other end the cut does not
//    split; chains of one curve check cut without a turn; chains along the
//    disc cuts, each with two offset variables behind the coordinates of the
//    copy points,
//  - the curve edges are chained on the cut copy (CurveChainBuilder).
class ChainBuilder {
public:
    struct Options {
        bool symmetry = false; // the part may be computed by halves
    };

    // Writes the result to part.chains (an earlier one is dropped). The part
    // has a cross field and a cut graph.
    static void compute(parts::Part& part, const Options& options, const run::Progress& progress);
};

} // namespace areko::remesher::param
