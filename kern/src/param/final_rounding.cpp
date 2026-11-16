#include "final_rounding.h"

namespace areko::remesher::param {

RoundingStep::Outcome FinalRounding::select(const SolveState& state, solver::DriverRequest& request)
{
    const solver::ConstrainedDriver& driver = state.driver;
    const int positions = static_cast<int>(driver.integerUnknowns().size());
    for (int q = 0; q < positions; ++q) {
        if (driver.isFree(q)) {
            request.selection.push_back({q, driver.nearestWhole(*request.x, q)});
        }
    }
    return Outcome::Delivered;
}

} // namespace areko::remesher::param
