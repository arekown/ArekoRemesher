#include "solution_builder.h"

#include "param/border_rounding.h"
#include "param/region_rounding.h"
#include "param/solution_seeds.h"
#include "param/solve_strategy.h"
#include "param/system_context.h"
#include "param/uv_output.h"

#include <utility>
#include <vector>

namespace areko::remesher::param {

void SolutionBuilder::compute(parts::Part& part, const Options& options, const run::Progress& progress)
{
    RegionRounding regions;
    BorderRounding borders;
    compute(part, options, progress, Steps{regions, borders});
}

void SolutionBuilder::compute(parts::Part& part, const Options& options, const run::Progress& progress, Steps steps)
{
    part.solution.reset();
    if (!part.system) {
        return;
    }
    System& system = *part.system;
    const SystemContext context(part, options.symmetry);

    SolveStrategy strategy(context, system, steps.regions, steps.borders, progress, options.report);
    solver::ConstrainedDriver driver(system.normal.matrix(), system.normal.rightHandSide(), strategy, progress);
    driver.setIntegerUnknowns(system.integerUnknowns, system.integerPeriods);
    strategy.attach(driver);
    std::vector<double> x;
    if (!driver.run(x)) {
        return;
    }
    progress.check();

    Solution solution;
    solution.rounds = strategy.rounds();
    solution.integerUnknowns = static_cast<int>(system.integerUnknowns.size());
    solution.uv = UvOutput::coordinates(system, x);
    prep::SymmetryCache* pairs = context.pairs();
    if (pairs != nullptr) {
        UvOutput::mirror(context, *pairs, solution.uv);
    }
    solution.singularities = part.chains->singularities;
    solution.seeds = part.chains->curves.turnPoints;
    SolutionSeeds::addBorderPoints(context, options.freezeBorders, solution);
    progress.check();
    if (pairs != nullptr) {
        SolutionSeeds::mirror(context, *pairs, solution);
    }
    if (options.report) {
        options.report(1.0f);
    }
    part.solution = std::move(solution);
}

} // namespace areko::remesher::param
