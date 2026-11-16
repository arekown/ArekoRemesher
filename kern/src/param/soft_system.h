#pragma once

#include "param/constraint_set.h"
#include "param/system.h"
#include "solver/sparse_matrix.h"

#include <span>
#include <vector>

namespace areko::remesher::run {
class Progress;
}

namespace areko::remesher::param {

// The additional system once its hard rows are eliminated: a small least
// squares system over its own unknowns.
struct ReducedSoftSystem {
    // The soft rows over the unknowns, then one row per unknown that ties it
    // weakly to its present value.
    solver::ListMatrix rows;
    std::vector<double> weights;
    std::vector<double> targets;
    // Column -> unknown (-1: the column is in no soft row) and back; the
    // unknowns are numbered as the soft rows name their columns.
    std::vector<int> unknownOf;
    std::vector<int> columnOf;
    // The columns the hard rows eliminated, each as a sum over other columns;
    // a term of the constant column stands for minus its coefficient.
    ConstraintSet conditions;

    [[nodiscard]] int unknownCount() const noexcept { return static_cast<int>(columnOf.size()); }
};

// The additional system of the rounding by regions. Its columns are the
// variables of the main system, one column for constants and one per bridge
// for the number of grid lines the bridge spans. Soft rows ask for
// "sum of coefficient * column = target" with a weight; hard rows have to
// hold exactly.
class SoftSystem {
public:
    SoftSystem(int variableCount, int bridgeCount)
        : variableCount_(variableCount), columnCount_(variableCount + 1 + bridgeCount)
    {
    }

    [[nodiscard]] int variableCount() const noexcept { return variableCount_; }
    [[nodiscard]] int constantColumn() const noexcept { return variableCount_; }
    [[nodiscard]] int firstBridgeColumn() const noexcept { return variableCount_ + 1; }
    [[nodiscard]] int bridgeColumn(int bridge) const noexcept { return firstBridgeColumn() + bridge; }
    [[nodiscard]] int columnCount() const noexcept { return columnCount_; }

    // A row counts even when it has no terms.
    void addSoft(TermList row, double target, double weight);
    void addHard(TermList row, double target);

    [[nodiscard]] const std::vector<TermList>& softRows() const noexcept { return soft_; }
    [[nodiscard]] const std::vector<double>& softTargets() const noexcept { return softTargets_; }
    [[nodiscard]] const std::vector<double>& softWeights() const noexcept { return softWeights_; }
    [[nodiscard]] const std::vector<TermList>& hardRows() const noexcept { return hard_; }
    [[nodiscard]] const std::vector<double>& hardTargets() const noexcept { return hardTargets_; }

    // Uses the rows up:
    //  - the conditions of the main system are put into all rows,
    //  - the hard rows are eliminated (their targets go to the constant
    //    column; `preferred` variables become pivots first) and the result is
    //    put into the soft rows,
    //  - the columns left in the soft rows are numbered in the order the
    //    rows name them,
    //  - every such column gets a row with a small weight: a variable of the
    //    main system is tied to its value in the solution `x` (0 without
    //    unknown), every other column to 1.
    // A stop request throws from the progress check.
    [[nodiscard]] ReducedSoftSystem reduce(const System& system, std::span<const int> preferred,
                                           std::span<const double> x, const run::Progress& progress) &&;

private:
    int variableCount_;
    int columnCount_;
    std::vector<TermList> soft_;
    std::vector<double> softTargets_;
    std::vector<double> softWeights_;
    std::vector<TermList> hard_;
    std::vector<double> hardTargets_;
};

} // namespace areko::remesher::param
