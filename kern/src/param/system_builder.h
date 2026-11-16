#pragma once

#include "param/elimination_columns.h"
#include "param/system.h"
#include "param/system_context.h"
#include "parts/part.h"
#include "run/progress.h"

#include <span>
#include <vector>

namespace areko::remesher::param {

// The set-up of the system of one part, step by step. `run` takes the steps
// in their order; a caller that brings rows of its own takes them one by one.
// The part has a cross field, a cut graph and chains and outlives the set-up,
// and so does the progress.
class SystemSetup {
public:
    // How a set-up ends.
    enum class Outcome {
        Built,              // the system stands
        EmptyMatrix,        // its normal matrix has no entry
        TermWithoutUnknown  // a row keeps a variable that is no unknown
    };

    // `symmetry`: the part may be computed by halves.
    SystemSetup(parts::Part& part, bool symmetry, const run::Progress& progress);

    // All steps. Only a set-up that ends with `Built` has a system to take.
    [[nodiscard]] Outcome run();

    // Variables with their kinds, the anchors and the transitions across the
    // chains: inside the chains as conditions, at their ends eliminated.
    void constrainChains();

    // Runs an elimination over `rows` and appends its conditions; they are
    // then put into all earlier conditions, with `intoAll` also into
    // themselves. Returns whether the elimination succeeded: `rows` is used
    // up then. After a failure the rows that are left and the conditions
    // found so far speak of variables like all others, so that the rows can
    // go on as rows of the system; a condition that was not reached is
    // appended empty and is the condition of no variable.
    bool eliminate(std::vector<TermList>& rows, EliminationMode mode, std::span<const int> preferred, bool intoAll);

    // Energy rows, behind them `penalties` (rows "= 0" over variables) with a
    // very large weight; all conditions put into the rows; the unknowns
    // numbered, the rows rewritten to them, the normal equations built and
    // the unknowns listed that have to become whole numbers.
    [[nodiscard]] Outcome assemble(std::vector<TermList> penalties);

    [[nodiscard]] const System& system() const noexcept { return system_; }
    [[nodiscard]] System take() noexcept;

private:
    // Iso lines along the curve chains. Returns the rows that could not be
    // eliminated.
    [[nodiscard]] std::vector<TermList> constrainIsoLines();

    SystemContext context_;
    const run::Progress& progress_;
    System system_;
    // First condition that came from a pivot other than 1: from there on the
    // grid periods have to be made to fit (-1: none).
    int periodStart_ = -1;
};

// Sets up the system of the parameterisation of a part from its chains, its
// combed cross field and its size field.
//  - variables: two coordinates per point of the cut copy and two offsets per
//    chain along a disc cut, each marked continuous or to become whole,
//  - conditions, in this order: the anchors; the transitions inside the
//    chains; the transitions at the chain ends, eliminated; the iso lines
//    along the curve chains, first those that tie a free variable to the
//    reference of its chain, then the others, eliminated. After each step
//    the new conditions are put into all earlier ones, so that in the end
//    every condition speaks of free variables only,
//  - the grid periods are made to fit the conditions,
//  - energy: four rows per triangle; rows of iso lines that could not be
//    eliminated follow them with a very large weight; all conditions are put
//    into the rows,
//  - the free variables are numbered as unknowns, the rows rewritten to
//    them, the normal equations built and the unknowns listed that have to
//    become whole numbers.
class SystemBuilder {
public:
    struct Options {
        bool symmetry = false; // the part may be computed by halves
    };

    // Writes the result to part.system (an earlier one is dropped). The part
    // has a cross field, a cut graph and chains. The system stays empty when
    // its normal matrix has no entry or a row keeps a variable that is no
    // unknown: the part cannot be solved then, and the answer says why.
    static SystemSetup::Outcome compute(parts::Part& part, const Options& options, const run::Progress& progress);
};

} // namespace areko::remesher::param
