#pragma once

#include "trace.h"

#include "parts/part.h"
#include "run/progress.h"

namespace areko::remesher::extract {

// Traces the whole-numbered grid lines of a part's parameterisation over its
// triangles and walks around every face of the grid:
//  - first from the start points (the singular points, else the other start
//    points of the solution): from every copy of a start point with whole
//    coordinates, in every triangle at it, one walk with u growing and one
//    with v growing,
//  - when that closes no face, the border points of the cut copy with whole
//    coordinates are added to the start points once and the first pass is
//    dropped and run again,
//  - then from every node made on the way, newest first, in each of the four
//    directions not done yet,
//  - finally the walks put back during that phase.
class IsolineTracer {
public:
    // Writes the result to part.trace (an earlier one is dropped). The part
    // stays without trace when it has no cut graph or solution, when no edge
    // has a remaining jump or when there is no start point. Every closed face
    // goes to `closer`. A stop request throws from the progress check.
    static void trace(parts::Part& part, FaceCloser& closer, const run::Progress& progress);
};

} // namespace areko::remesher::extract
