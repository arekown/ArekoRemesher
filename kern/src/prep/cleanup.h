#pragma once

#include "preparation.h"

#include "run/progress.h"

namespace areko::remesher::prep {

// Second part of the preparation (stage "cleaned"), on the result of the first:
// denoising, cleaning, sliver and small region treatment with the line care
// around it, the second decimation to the work resolution, a second cleaning
// and the refinement of coarse meshes.
class Cleanup {
public:
    static void run(Prepared& prepared, run::Progress& progress);
};

} // namespace areko::remesher::prep
