#include "elimination.h"

#include "param/constraint_terms.h"
#include "run/progress.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace areko::remesher::param {

namespace {

// Columns handled between two cancel checks.
constexpr int blockSize = 1000;
// A pivot of smaller magnitude cannot be divided by.
constexpr double pivotLimit = static_cast<double>(1e-10f);
static_assert(pivotLimit == 1.000000013351432e-10);
// Terms below this magnitude are dropped from the conditions of the soft system.
constexpr double tinyLimit = 1e-8;

[[nodiscard]] bool isUnit(double value) noexcept
{
    return std::fabs(value) == 1.0;
}

[[nodiscard]] bool leadsWith(const TermList& row, int column) noexcept
{
    return !row.empty() && row.front().index == column;
}

// Position of the first entry that fulfils `wanted`, looking only at the
// entries in front of column `limit`; -1 if there is none.
template <class Wanted>
[[nodiscard]] int firstEntryBelow(const TermList& row, int limit, Wanted wanted)
{
    for (std::size_t i = 0; i < row.size(); ++i) {
        if (row[i].index >= limit) {
            break;
        }
        if (wanted(row[i].value)) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

} // namespace

Elimination::Elimination(std::vector<TermList>& rows, int variableCount, const EliminationSetup& setup,
                         const run::Progress& progress)
    : rows_(rows),
      progress_(progress),
      mode_(setup.mode),
      flags_(setup.flags),
      rowCount_(static_cast<int>(rows.size())),
      columns_(rows, variableCount, setup)
{
    for (TermList& row : rows_) {
        solver::ListMatrix::order(row, solver::ListMatrix::SameIndex::Separate);
    }
    continuousEnd_ = continuousEnd();
    columnToConstraint_.assign(static_cast<std::size_t>(columns_.count()), -1);
}

// Number of leading columns whose variables are continuous (all columns, if
// every variable is). Up to here a continuous column may be exchanged.
int Elimination::continuousEnd() const
{
    if (mode_ == EliminationMode::SoftSystem) {
        return 0;
    }
    int end = 0;
    while (end < columns_.count()
           && flags_[static_cast<std::size_t>(columns_.variableOf(end))] == EliminationColumns::continuousFlag) {
        ++end;
    }
    return end;
}

bool Elimination::run()
{
    if (!ran_) {
        ran_ = true;
        progress_.check();
        forward();
        progress_.check();
        eliminated_ = substitute();
    }
    return eliminated_;
}

void Elimination::forward()
{
    const bool soft = mode_ == EliminationMode::SoftSystem;
    // The column of the constants is the last one and never a pivot column.
    const int columnEnd = columns_.count() - (soft ? 1 : 0);
    if (columnEnd < 1) {
        return;
    }
    int row = 0;
    int column = 0;
    int blockEnd = std::min(blockSize, columnEnd);
    while (true) {
        progress_.check();
        while (column < blockEnd) {
            if (soft && row >= rowCount_) {
                break;
            }
            if (soft ? pickInRow(row, column) : pickInColumn(row, column)) {
                // With a pivot in every row only the current block of columns
                // ends. A further block goes on with the same row and column:
                // the last pivot row is taken once more and the conditions
                // end with a copy of the last one.
                if (pivotCount() == rowCount_) {
                    break;
                }
                ++row;
                ++column;
            } else if (soft) {
                ++row;
            } else {
                ++column;
            }
        }
        if (blockEnd >= columnEnd) {
            break;
        }
        blockEnd = std::min(blockEnd + blockSize, columnEnd);
    }
}

// Pivot for `column` among the rows from `row` on: a row that leads with a
// coefficient of magnitude 1 in this column; else a row with a coefficient of
// magnitude 1 in an exchangeable column, which is exchanged with this one;
// else the row that leads in this column with the smallest magnitude.
bool Elimination::pickInColumn(int row, int column)
{
    int found = -1;
    for (int s = row; s < rowCount_; ++s) {
        if (leadsWith(rowAt(s), column) && isUnit(rowAt(s).front().value)) {
            found = s;
            break;
        }
    }
    if (found < 0) {
        const bool continuous = flags_[static_cast<std::size_t>(columns_.variableOf(column))]
                                == EliminationColumns::continuousFlag;
        const int limit = continuous ? continuousEnd_ : columns_.split();
        for (int s = row; s < rowCount_; ++s) {
            const int entry = firstEntryBelow(rowAt(s), limit, isUnit);
            if (entry >= 0) {
                columns_.swap(rows_, column, rowAt(s)[static_cast<std::size_t>(entry)].index);
                found = s;
                break;
            }
        }
    }
    if (found < 0) {
        // Strictly smaller wins, so the first row on ties; 0 and NaN never.
        double best = ConstraintTerms::magnitudeStart;
        for (int s = row; s < rowCount_; ++s) {
            if (!leadsWith(rowAt(s), column)) {
                continue;
            }
            const double magnitude = std::fabs(rowAt(s).front().value);
            if (magnitude != 0.0 && magnitude < best) {
                best = magnitude;
                found = s;
            }
        }
        if (found >= 0 && firstSmallestPivot_ == -1) {
            firstSmallestPivot_ = row;
        }
    }
    if (found < 0) {
        return false;
    }
    pivot(found, row, column);
    return true;
}

// Pivot for `column` in the given row: its leading coefficient of magnitude 1
// in this column; else its first coefficient of magnitude 1, else its first
// coefficient that is not 0, in front of the split, whose column is exchanged
// with this one.
bool Elimination::pickInRow(int row, int column)
{
    const TermList& entries = rowAt(row);
    if (!(leadsWith(entries, column) && isUnit(entries.front().value))) {
        const int split = columns_.split();
        int entry = firstEntryBelow(entries, split, isUnit);
        if (entry < 0) {
            entry = firstEntryBelow(entries, split, [](double value) { return value != 0.0; });
        }
        if (entry < 0) {
            return false;
        }
        columns_.swap(rows_, column, entries[static_cast<std::size_t>(entry)].index);
    }
    pivot(row, row, column);
    return true;
}

// Moves the pivot row to `row` and removes the column from every later row
// that leads with it: later row * pivot - pivot row * its coefficient.
void Elimination::pivot(int found, int row, int column)
{
    if (found != row) {
        std::swap(rowAt(found), rowAt(row));
    }
    const TermList& pivotRow = rowAt(row);
    const double pivotValue = pivotRow.front().value;
    for (int s = row + 1; s < rowCount_; ++s) {
        TermList& later = rowAt(s);
        if (leadsWith(later, column) && later.front().value != 0.0) {
            later = ConstraintTerms::combine(-later.front().value, pivotRow, pivotValue, later);
            ConstraintTerms::reduceToIntegers(later);
        }
    }
    columnToConstraint_[static_cast<std::size_t>(column)] = pivotCount();
    pivotColumns_.push_back(column);
    pivotRows_.push_back(row);
}

// From the last pivot backwards: the pivot row divided by its pivot, with the
// conditions found so far put in for their columns.
bool Elimination::substitute()
{
    const bool soft = mode_ == EliminationMode::SoftSystem;
    constraints_.assign(pivotColumns_.size(), Constraint{});
    bool dropped = false;
    TermList terms;
    for (std::size_t k = pivotColumns_.size(); k-- > 0;) {
        const TermList& row = rowAt(pivotRows_[k]);
        if (row.empty()) {
            continue;
        }
        const double pivotValue = row.front().value;
        if (pivotLimit > std::fabs(pivotValue)) {
            if (!soft) {
                return false;
            }
            dropped = true;
            continue;
        }
        terms.clear();
        for (std::size_t j = 1; j < row.size(); ++j) {
            const double factor = -row[j].value / pivotValue;
            const int source = columnToConstraint_[static_cast<std::size_t>(row[j].index)];
            if (source == -1) {
                terms.push_back({row[j].index, factor});
            } else {
                for (const Term& term : constraints_[static_cast<std::size_t>(source)].terms) {
                    terms.push_back({term.index, factor * term.value});
                }
            }
        }
        ConstraintTerms::normalize(terms);
        if (soft) {
            ConstraintTerms::dropTiny(terms, tinyLimit);
        }
        constraints_[k].terms = terms;
        constraints_[k].pivot = pivotColumns_[k];
    }
    rows_.clear();
    progress_.check();
    restoreVariables();
    if (dropped) {
        std::erase_if(constraints_, [](const Constraint& constraint) { return constraint.pivot < 0; });
    }
    return true;
}

// Column -> variable for every term and every pivot; the terms keep the order
// of their columns.
void Elimination::restoreVariables()
{
    if (!ran_ || !inColumns_) {
        return;
    }
    inColumns_ = false;
    const auto restore = [&](TermList& terms) {
        for (Term& term : terms) {
            term.index = columns_.variableOf(term.index);
        }
    };
    for (TermList& row : rows_) {
        restore(row);
    }
    for (Constraint& constraint : constraints_) {
        if (constraint.pivot >= 0) {
            constraint.pivot = columns_.variableOf(constraint.pivot);
        }
        restore(constraint.terms);
    }
}

void Elimination::enterPivots(std::vector<int>& variableToConstraint, int firstIndex) const
{
    for (std::size_t i = 0; i < constraints_.size(); ++i) {
        const int pivot = constraints_[i].pivot;
        if (pivot >= 0) {
            variableToConstraint[static_cast<std::size_t>(pivot)] = firstIndex + static_cast<int>(i);
        }
    }
}

} // namespace areko::remesher::param
