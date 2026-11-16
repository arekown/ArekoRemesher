#pragma once

#include "preparation.h"

#include "run/progress.h"

namespace areko::remesher::prep {

// Symmetry part of the preparation, after the refinement:
//  - only with symmetry on: lines into the sharp edge attribute (when the
//    input had one), then one cut per active axis X, Y, Z,
//  - one axis and no guiding curve that cuts into the mesh: cut, mirror and
//    keep the mirror pairs; otherwise every axis only cuts, the half mode is
//    set and the quad budget becomes the quad count / 2^axes,
//  - the plane position of every cut axis is kept; lines again from the
//    sharp edges,
//  - a polygon mesh is triangulated (mirror path when an axis was cut) and
//    repaired when edges with more than two corners remain; the mirror pairs
//    are kept as they are.
class Symmetry {
public:
    static void run(Prepared& prepared, run::Progress& progress);
};

} // namespace areko::remesher::prep
