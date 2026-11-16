#include "constraint_periods.h"

#include "util/rounding.h"

#include <cmath>
#include <numeric>

namespace areko::remesher::param {

namespace {

constexpr double wholeTolerance = 1e-7;
constexpr int factorEnd = 50;

// Whole within the tolerance: a whole number lies between x - tolerance and
// x + tolerance, or x is 0 within the tolerance. NaN, infinite values and
// values outside the int range are never whole.
[[nodiscard]] bool isWhole(double x) noexcept
{
    return util::truncateToInt(x - wholeTolerance) != util::truncateToInt(x + wholeTolerance)
           || std::fabs(x) <= wholeTolerance;
}

} // namespace

int ConstraintPeriods::adjust(std::vector<int>& periods, std::span<const Constraint> constraints, int first)
{
    int changed = 0;
    if (first < 0) {
        return changed;
    }
    for (std::size_t b = static_cast<std::size_t>(first); b < constraints.size(); ++b) {
        const Constraint& constraint = constraints[b];
        if (constraint.pivot < 0) {
            continue;
        }
        const int pivotPeriod = periods[static_cast<std::size_t>(constraint.pivot)];
        if (pivotPeriod == 0) {
            continue;
        }
        for (const Term& term : constraint.terms) {
            int& period = periods[static_cast<std::size_t>(term.index)];
            const int old = period;
            const double share = static_cast<double>(old) * term.value / static_cast<double>(pivotPeriod);
            if (isWhole(share)) {
                continue;
            }
            int factor = 2;
            while (factor < factorEnd && !isWhole(static_cast<double>(factor) * share)) {
                ++factor;
            }
            if (factor == factorEnd) {
                continue;
            }
            const auto wide = static_cast<unsigned>(factor);
            const auto before = static_cast<unsigned>(old);
            period = static_cast<int>(wide * before / std::gcd(wide, before));
            if (old == 1) {
                ++changed;
            }
        }
    }
    return changed;
}

} // namespace areko::remesher::param
