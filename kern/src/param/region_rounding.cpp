#include "region_rounding.h"

#include "param/bridge_rows.h"
#include "param/region_bridges.h"
#include "param/region_growth.h"
#include "param/seam_rows.h"
#include "param/soft_system.h"
#include "param/uv_output.h"
#include "solver/normal_equations.h"
#include "util/rounding.h"

#include <array>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

namespace areko::remesher::param {

using solver::DriverRequest;

namespace {

// A soft row per offset across a disc cut that is an unknown of the main
// system: the offset keeps its present value.
void addOffsetRows(const SolveState& state, std::span<const double> x, SoftSystem& system)
{
    const Chains& chains = state.context.chains();
    for (const int chain : chains.genusChains) {
        const int first = chains.offsetVariable[static_cast<std::size_t>(chain)];
        for (const int variable : {first, first + 1}) {
            const int unknown = state.system.variableToUnknown[static_cast<std::size_t>(variable)];
            if (unknown != -1) {
                system.addSoft(TermList{{variable, 1.0}}, x[static_cast<std::size_t>(unknown)], 1.0);
            }
        }
    }
}

// The variables the elimination of the hard rows turns into pivots first:
// with grid periods in the main solve those of its integer list whose period
// is 1.
[[nodiscard]] std::vector<int> preferredPivots(const SolveState& state)
{
    const std::vector<int>& unknowns = state.driver.integerUnknowns();
    const std::vector<double>& periods = state.driver.periods();
    std::vector<int> preferred;
    for (std::size_t k = 0; k < periods.size(); ++k) {
        if (periods[k] == 1.0) {
            const int variable = state.system.unknownToVariable[static_cast<std::size_t>(unknowns[k])];
            if (variable >= 0) {
                preferred.push_back(variable);
            }
        }
    }
    return preferred;
}

// The unknowns of the additional system that have to become whole numbers,
// with their periods: the variables of the integer list of the main solve,
// in its order, then every column behind the variables.
struct WholeUnknowns {
    std::vector<int> unknowns;
    std::vector<double> periods;
};

[[nodiscard]] WholeUnknowns wholeUnknowns(const SolveState& state, const ReducedSoftSystem& reduced)
{
    WholeUnknowns out;
    for (const int mainUnknown : state.driver.integerUnknowns()) {
        const auto variable = static_cast<std::size_t>(state.system.unknownToVariable[static_cast<std::size_t>(mainUnknown)]);
        if (reduced.unknownOf[variable] != -1) {
            out.unknowns.push_back(reduced.unknownOf[variable]);
            out.periods.push_back(static_cast<double>(state.system.periods[variable]));
        }
    }
    for (std::size_t column = static_cast<std::size_t>(state.system.variableCount); column < reduced.unknownOf.size();
         ++column) {
        if (reduced.unknownOf[column] >= 0) {
            out.unknowns.push_back(reduced.unknownOf[column]);
            out.periods.push_back(1.0);
        }
    }
    return out;
}

// Rounds of the solve of the additional system: in the first round every
// bridge count is rounded at once, to at least 1; everything else is left to
// the batch rounding of the driver.
class BridgeCounts final : public solver::DriverCallback {
public:
    // `unknowns` is the integer list of the solve, `firstBridgeColumn` the
    // first column of a bridge count.
    BridgeCounts(std::span<const int> unknowns, std::span<const int> columnOf, int firstBridgeColumn)
        : unknowns_(unknowns), columnOf_(columnOf), firstBridgeColumn_(firstBridgeColumn)
    {
    }

    void onRequest(DriverRequest& request) override
    {
        if (asked_) {
            return;
        }
        asked_ = true;
        const std::vector<double>& y = *request.x;
        for (std::size_t k = 0; k < unknowns_.size(); ++k) {
            const auto unknown = static_cast<std::size_t>(unknowns_[k]);
            if (columnOf_[unknown] >= firstBridgeColumn_) {
                const int count = util::roundHalfAway(y[unknown]);
                request.selection.push_back({static_cast<int>(k), count > 0 ? static_cast<double>(count) : 1.0});
            }
        }
        if (!request.selection.empty()) {
            request.answer = DriverRequest::Answer::Delivered;
        }
    }

private:
    std::span<const int> unknowns_;
    std::span<const int> columnOf_;
    int firstBridgeColumn_;
    bool asked_ = false;
};

// Hands the solved additional system `y` back to the integer list of the
// main solve: a variable that is one of its unknowns takes its value; a
// variable one of its conditions eliminated takes the value of the
// condition, rounded to its period. Returns whether a position it has a
// value for was left out because it is fixed already.
[[nodiscard]] bool transferBack(const SolveState& state, const ReducedSoftSystem& reduced, int constantColumn,
                                std::span<const double> y, DriverRequest& request)
{
    const std::vector<int>& unknowns = state.driver.integerUnknowns();
    const std::vector<double>& periods = state.driver.periods();
    const mesh::BitArray& fixed = state.driver.fixedPositions();
    const std::vector<Constraint>& conditions = reduced.conditions.list();
    // Position of the integer list every condition stands for, -1: none.
    std::vector<int> positionOf(conditions.size(), -1);
    bool leftOut = false;
    for (std::size_t k = 0; k < unknowns.size(); ++k) {
        const int variable = state.system.unknownToVariable[static_cast<std::size_t>(unknowns[k])];
        const int unknown = reduced.unknownOf[static_cast<std::size_t>(variable)];
        const int condition = unknown == -1 ? reduced.conditions.ofVariable(variable) : -1;
        if (unknown == -1 && condition == -1) {
            continue;
        }
        // A position that is fixed already keeps its value.
        if (fixed.testBounded(k)) {
            leftOut = true;
            continue;
        }
        if (unknown != -1) {
            request.selection.push_back({static_cast<int>(k), y[static_cast<std::size_t>(unknown)]});
        } else {
            positionOf[static_cast<std::size_t>(condition)] = static_cast<int>(k);
        }
    }
    for (std::size_t b = conditions.size(); b-- > 0;) {
        double value = 0.0;
        for (const Term& term : conditions[b].terms) {
            if (term.index == constantColumn) {
                value = value - term.value;
            } else if (const int unknown = reduced.unknownOf[static_cast<std::size_t>(term.index)]; unknown >= 0) {
                value = value + term.value * y[static_cast<std::size_t>(unknown)];
            }
        }
        const int position = positionOf[b];
        if (position != -1) {
            const double period = periods.empty() ? 1.0 : periods[static_cast<std::size_t>(position)];
            request.selection.push_back({position, solver::ConstrainedDriver::toPeriod(value, period).value});
        }
    }
    return leftOut;
}

} // namespace

RoundingStep::Outcome RegionRounding::select(const SolveState& state, DriverRequest& request)
{
    const SystemContext& context = state.context;
    const System& system = state.system;
    const std::vector<double>& x = *request.x;
    state.progress.check();

    const UvOutput::Coordinates uv = UvOutput::coordinates(system, x);
    const std::array<RegionGrowth, 2> growth{RegionGrowth(context, uv, 0), RegionGrowth(context, uv, 1)};
    RegionBridges bridges;
    for (const RegionGrowth& regions : growth) {
        bridges.add(regions, uv, context.copy().adjacency());
    }
    state.progress.check();

    SoftSystem rows(system.variableCount, bridges.count());
    BridgeRows::add(context, system.constraints, uv, growth, bridges, rows);
    addOffsetRows(state, x, rows);
    SeamRows::add(context, uv, rows);
    state.progress.check();

    const int constantColumn = rows.constantColumn();
    const int firstBridgeColumn = rows.firstBridgeColumn();
    const std::vector<int> preferred = preferredPivots(state);
    const ReducedSoftSystem reduced = std::move(rows).reduce(system, preferred, x, state.progress);
    WholeUnknowns whole = wholeUnknowns(state, reduced);
    solver::NormalEquations::ListSystem normal = solver::NormalEquations::assembleLists(
        reduced.rows, reduced.weights, reduced.targets, state.progress);
    state.progress.check();

    BridgeCounts counts(whole.unknowns, reduced.columnOf, firstBridgeColumn);
    solver::DriverOptions options;
    options.batchMode = true;
    options.batchBudget = batchBudget;
    solver::ConstrainedDriver driver(normal.matrix, normal.rhs, counts, state.progress, options);
    // The callback reads the list through its span; the driver gets a copy.
    driver.setIntegerUnknowns(whole.unknowns, std::move(whole.periods));
    std::vector<double> y;
    if (!driver.run(y)) {
        return Outcome::Failed;
    }
    const bool leftOut = transferBack(state, reduced, constantColumn, y, request);
    return leftOut && request.selection.empty() ? Outcome::AlreadyFixed : Outcome::Delivered;
}

} // namespace areko::remesher::param
