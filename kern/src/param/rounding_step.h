#pragma once

#include "param/system.h"
#include "param/system_context.h"
#include "run/progress.h"
#include "solver/constrained_driver.h"

namespace areko::remesher::param {

// What a rounding step reads while the system of a part is being solved. The
// latest solution is request.x of the round.
struct SolveState {
    const SystemContext& context;            // part, mesh, cut copy, cut graph, chains, symmetry
    const System& system;                    // kinds, periods, conditions, numbering, rows
    const solver::ConstrainedDriver& driver; // integer list, periods and fixed positions of the solve
    const run::Progress& progress;
};

// A step of the solving phase that proposes values for unknowns that have to
// become whole numbers: it appends (position in the integer list, value) to
// the selection of the round; the driver fixes them and solves again. A step
// changes neither the system nor the driver.
class RoundingStep {
public:
    enum class Outcome {
        Failed,      // the step found no values
        Delivered,   // the selection holds what the step found
        AlreadyFixed // it found values, but only for positions fixed already
    };

    virtual ~RoundingStep() = default;
    virtual Outcome select(const SolveState& state, solver::DriverRequest& request) = 0;
};

} // namespace areko::remesher::param
