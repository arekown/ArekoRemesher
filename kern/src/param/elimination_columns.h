#pragma once

#include "param/constraint.h"

#include <cstdint>
#include <span>
#include <vector>

namespace areko::remesher::param {

// What the rows of an elimination are and how their columns are ordered.
enum class EliminationMode {
    ChainEnds,   // transitions at the ends of the cut chains
    IsoLines,    // iso lines along curves: continuous variables come first
    SoftSystem   // hard rows of the additional system, with a column of constants
};

struct EliminationSetup {
    EliminationMode mode = EliminationMode::ChainEnds;
    std::span<const int> preferred;        // variables that get their columns before those of the rows
    std::span<const std::uint8_t> flags;   // one per variable (ChainEnds, IsoLines)
    int constantColumn = -1;               // SoftSystem: variable that carries the constants
};

// Column order of an elimination: every variable that occurs gets a column,
// in the order in which it should become a pivot. Columns from `split` on are
// never swapped into a pivot position.
class EliminationColumns {
public:
    static constexpr std::uint8_t continuousFlag = 0;
    static constexpr std::uint8_t keptFreeFlag = 2;   // these variables get the last columns

    // Numbers the variables and rewrites every entry of the rows to its column:
    //  - IsoLines only: variables with the continuous flag, in order of first occurrence,
    //  - the preferred variables in their order (also those that occur in no row),
    //  - ChainEnds, IsoLines: first occurrences of all variables without the
    //    kept-free flag (up to here: split), then the remaining ones,
    //  - SoftSystem: variables below the constant column, then those above it
    //    (up to here: split), last the constant column itself.
    EliminationColumns(std::vector<TermList>& rows, int variableCount, const EliminationSetup& setup);

    [[nodiscard]] int count() const noexcept { return static_cast<int>(variableOf_.size()); }
    [[nodiscard]] int split() const noexcept { return split_; }
    [[nodiscard]] int variableOf(int column) const { return variableOf_[static_cast<std::size_t>(column)]; }
    [[nodiscard]] std::span<const int> variables() const noexcept { return variableOf_; }

    // Exchanges two columns: in every row, which is sorted again, and in the numbering.
    void swap(std::vector<TermList>& rows, int a, int b);

private:
    std::vector<int> variableOf_;
    int split_ = 0;
};

} // namespace areko::remesher::param
