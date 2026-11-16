#pragma once

#include "run/progress.h"

namespace areko::remesher::parts {
struct Part;
} // namespace areko::remesher::parts

namespace areko::remesher::extract {

// The extraction of one part, in two stages:
//  - raw mesh: the grid lines are traced and every closed face becomes
//    polygons, the raw mesh is built from them, its border classified, then
//    degenerate faces are cleaned, doublets dissolved and holes, border loops
//    and pinched points closed,
//  - quads: the faces of five or more corners of the raw mesh are cut into
//    quads; the point data grow with the new points (no source, no triangle,
//    no neighbours).
class Extraction {
public:
    struct Options {
        bool followBorders = true;
        // Carry the point colours of the triangle mesh over to the result.
        bool keepColors = false;
        // The symmetry is switched on: the edges on the mirror plane of a
        // symmetric part stay fixed in the final relaxation.
        bool symmetry = false;
        // Crease chains are cut into the quad mesh (otherwise only the edges
        // between their points become curve edges).
        bool cutAlongCurves = true;
    };

    // Writes the raw mesh to part.raw (an earlier one is dropped); the part
    // stays without when the trace finds nothing to trace. A stop request
    // throws from the progress check.
    static void rawMesh(parts::Part& part, const Options& options, const run::Progress& progress);

    // Turns part.raw into the quad mesh (nothing without raw mesh).
    static void quads(parts::Part& part, const Options& options);

    // Brings crease curves, sharp edges and face groups of the triangle mesh
    // onto the quad mesh (after quads(); nothing without raw mesh).
    static void features(parts::Part& part, const Options& options);

    // Finishes the quad mesh of part.raw (nothing without raw mesh), after
    // the features are carried over: the points are relaxed and walked back
    // onto the triangle mesh, the kept border points are put back exactly
    // with frozen borders, the point colours are carried over when asked
    // for, and degenerate faces and unused points are cleaned. A stop
    // request throws from the progress check.
    static void finish(parts::Part& part, const Options& options, const run::Progress& progress);
};

} // namespace areko::remesher::extract
