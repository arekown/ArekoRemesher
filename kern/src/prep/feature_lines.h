#pragma once

#include "preparation.h"

#include "run/progress.h"

namespace areko::remesher::prep {

// Feature lines, the last part of the preparation (after the symmetry):
//  - final lines from the sharp edges and the face group borders, and the
//    point map of their points,
//  - with denoise or crease detection on: crease chains on the mesh and the
//    smoothing after the detection; the area is stored when none was,
//  - with crease detection on: face group lines that stand for a crease are
//    marked, then the detected chains join the lines (type = class, flag 4);
//    without it the detected chains are dropped.
class FeatureLines {
public:
    static void run(Prepared& prepared, run::Progress& progress);
};

} // namespace areko::remesher::prep
