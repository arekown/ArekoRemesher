#pragma once

#include "param/constraint.h"
#include "param/elimination_columns.h"

#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace areko::remesher::run {
class Progress;
}

namespace areko::remesher::param {

// Eliminates linear conditions "row = 0" without fractions: columns are taken
// in a fixed order, a coefficient of magnitude 1 is preferred as pivot (if
// necessary by exchanging columns), rows are combined with whole factors, and
// the pivots are then expressed by the remaining variables from the last one
// backwards. The result is one condition per pivot, x[pivot] = sum over
// variables that are not eliminated.
//
// The object prepares the rows on construction (their entries are renumbered
// to columns and sorted) and is run once. The rows, the flags of the setup
// and the progress must outlive it; the preferred variables are read on
// construction only.
class Elimination {
public:
    Elimination(std::vector<TermList>& rows, int variableCount, const EliminationSetup& setup,
                const run::Progress& progress);

    // True: the rows are used up (the list is empty) and the conditions carry
    // variable numbers. False (ChainEnds, IsoLines only; SoftSystem always
    // gives true): a pivot was too small; the rows stay as far as they were
    // eliminated, rows and the conditions found so far keep their column
    // numbers, the others stay empty. Cancelling throws from the progress
    // check. Only the first call computes: a further one gives the same
    // answer, false after a cancelled run.
    bool run();

    // After a run that gave false: writes the variables back into the rows
    // and into the conditions found so far, for a caller that goes on with
    // them. The terms keep the order of their columns. Does nothing where
    // the variables already stand.
    void restoreVariables();

    // The eliminated variables in the order of their pivots.
    [[nodiscard]] const std::vector<Constraint>& constraints() const noexcept { return constraints_; }
    [[nodiscard]] std::vector<Constraint> takeConstraints() noexcept { return std::move(constraints_); }

    // Number of the first condition whose pivot came from the rule of the
    // smallest magnitude (its coefficient may not be 1), or -1.
    [[nodiscard]] int firstSmallestPivot() const noexcept { return firstSmallestPivot_; }

    // Enters the conditions into a map variable -> condition as if they were
    // appended to a list at `firstIndex`. Conditions without a pivot get no entry.
    void enterPivots(std::vector<int>& variableToConstraint, int firstIndex) const;

private:
    [[nodiscard]] int continuousEnd() const;
    [[nodiscard]] int pivotCount() const noexcept { return static_cast<int>(pivotColumns_.size()); }
    [[nodiscard]] TermList& rowAt(int row) { return rows_[static_cast<std::size_t>(row)]; }
    void forward();
    bool pickInColumn(int row, int column);
    bool pickInRow(int row, int column);
    void pivot(int found, int row, int column);
    bool substitute();

    std::vector<TermList>& rows_;
    const run::Progress& progress_;
    EliminationMode mode_;
    std::span<const std::uint8_t> flags_;
    int rowCount_;
    EliminationColumns columns_;
    int continuousEnd_ = 0;

    std::vector<int> pivotColumns_;
    std::vector<int> pivotRows_;
    std::vector<int> columnToConstraint_;
    std::vector<Constraint> constraints_;
    int firstSmallestPivot_ = -1;
    bool ran_ = false;
    bool eliminated_ = false;
    // Rows and conditions carry column numbers.
    bool inColumns_ = true;
};

} // namespace areko::remesher::param
