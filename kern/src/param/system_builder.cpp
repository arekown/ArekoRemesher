#include "system_builder.h"

#include "param/anchor_constraints.h"
#include "param/constraint_periods.h"
#include "param/elimination.h"
#include "param/energy_rows.h"
#include "param/iso_line_constraints.h"
#include "param/transition_constraints.h"
#include "param/unknown_numbering.h"
#include "param/variable_flags.h"

#include <utility>

namespace areko::remesher::param {

namespace {

// Weight of the rows of iso lines that could not be eliminated.
constexpr double penaltyWeight = 1e20;

} // namespace

SystemSetup::SystemSetup(parts::Part& part, bool symmetry, const run::Progress& progress)
    : context_(part, symmetry), progress_(progress)
{
}

System SystemSetup::take() noexcept
{
    return std::move(system_);
}

SystemSetup::Outcome SystemSetup::run()
{
    constrainChains();
    std::vector<TermList> penalties = constrainIsoLines();
    system_.periodChanges = ConstraintPeriods::adjust(system_.periods, system_.constraints.list(), periodStart_);
    progress_.check();
    return assemble(std::move(penalties));
}

void SystemSetup::constrainChains()
{
    VariableFlags::assign(context_, system_);
    progress_.check();
    AnchorConstraints::add(context_, system_);
    TransitionConstraints transitions(context_, system_);
    transitions.addInterior();
    progress_.check();
    std::vector<TermList> rows = transitions.endRows();
    const std::vector<int> preferred = VariableFlags::markChainEnds(context_, system_);
    system_.chainEndsEliminated = eliminate(rows, EliminationMode::ChainEnds, preferred, true);
    system_.sections.chainEnds = system_.constraints.count();
    progress_.check();
}

bool SystemSetup::eliminate(std::vector<TermList>& rows, EliminationMode mode, std::span<const int> preferred,
                            bool intoAll)
{
    ConstraintSet& set = system_.constraints;
    EliminationSetup setup;
    setup.mode = mode;
    setup.preferred = preferred;
    setup.flags = system_.flags;
    Elimination elimination(rows, system_.variableCount, setup, progress_);
    const bool done = elimination.run();
    if (!done) {
        // What a failed run leaves is in the columns of the elimination.
        elimination.restoreVariables();
    }
    const int first = set.count();
    if (periodStart_ == -1 && elimination.firstSmallestPivot() != -1) {
        periodStart_ = first + elimination.firstSmallestPivot();
    }
    elimination.enterPivots(set.variableMap(), first);
    set.append(elimination.takeConstraints());
    set.resolve(first, intoAll ? set.count() : first);
    return done;
}

std::vector<TermList> SystemSetup::constrainIsoLines()
{
    ConstraintSet& set = system_.constraints;
    const int first = set.count();
    system_.sections.isoLines = first;
    if (context_.graph().curveEdges.empty()) {
        return {};
    }
    IsoLineConstraints::Result iso = IsoLineConstraints::build(context_, system_);
    progress_.check();
    set.resolve(first, first);
    progress_.check();
    set.substituteInto(iso.rows);
    system_.sections.isoLines = set.count();
    progress_.check();
    system_.isoLinesEliminated = eliminate(iso.rows, EliminationMode::IsoLines, iso.preferred, false);
    return std::move(iso.rows);
}

SystemSetup::Outcome SystemSetup::assemble(std::vector<TermList> penalties)
{
    EnergyRows::Result energy = EnergyRows::build(context_);
    progress_.check();
    system_.triangleRowStart = std::move(energy.triangleRowStart);
    system_.targets = std::move(energy.targets);
    system_.weights = std::move(energy.weights);
    system_.energyRowCount = static_cast<int>(energy.rows.size());
    system_.penaltyRowCount = static_cast<int>(penalties.size());
    // Rows that could not be eliminated are rows of the system like the
    // others: they are counted, renumbered and enter the normal equations.
    std::vector<TermList> rows = std::move(energy.rows);
    for (TermList& penalty : penalties) {
        rows.push_back(std::move(penalty));
        system_.targets.push_back(0.0);
        system_.weights.push_back(penaltyWeight);
    }
    system_.constraints.substituteInto(rows);
    progress_.check();

    UnknownNumbering::assign(context_, system_);
    system_.unresolvedTerms = UnknownNumbering::renumber(system_, rows);
    if (system_.unresolvedTerms > 0) {
        return Outcome::TermWithoutUnknown;
    }
    system_.rows = solver::ListMatrix(static_cast<int>(rows.size()), system_.unknownCount());
    for (std::size_t i = 0; i < rows.size(); ++i) {
        system_.rows.list(static_cast<int>(i)) = std::move(rows[i]);
    }
    if (system_.normal.assemble(system_.rows, system_.weights, system_.targets, solver::NormalEquations::Scope::All,
                                progress_) == solver::NormalEquations::Result::EmptyMatrix) {
        return Outcome::EmptyMatrix;
    }
    progress_.check();
    UnknownNumbering::listIntegers(system_);
    return Outcome::Built;
}

SystemSetup::Outcome SystemBuilder::compute(parts::Part& part, const Options& options, const run::Progress& progress)
{
    part.system.reset();
    SystemSetup setup(part, options.symmetry, progress);
    const SystemSetup::Outcome outcome = setup.run();
    if (outcome == SystemSetup::Outcome::Built) {
        part.system = setup.take();
    }
    return outcome;
}

} // namespace areko::remesher::param
