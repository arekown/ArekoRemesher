#pragma once

#include "param/constraint_set.h"
#include "param/elimination_columns.h"
#include "solver/normal_equations.h"
#include "solver/sparse_matrix.h"

#include <cstdint>
#include <vector>

namespace areko::remesher::param {

// What a variable of the system is, one value per variable. Every value but
// `continuous` marks a variable that has to become a whole number.
namespace variable_flag {
inline constexpr std::uint8_t continuous = EliminationColumns::continuousFlag;
inline constexpr std::uint8_t node = 1;          // copy of a singular or border node, point of a crease or border curve
inline constexpr std::uint8_t keptFree = EliminationColumns::keptFreeFlag; // at the end of a chain along a disc cut
inline constexpr std::uint8_t offset = 3;        // offset across a chain along a disc cut
inline constexpr std::uint8_t checkSource = 4;   // copy of a source of the curve check
inline constexpr std::uint8_t nodeCopy = 5;      // further copy of a node
inline constexpr std::uint8_t frozenBorder = 6;  // copy of a border point of a frozen border
} // namespace variable_flag

// Where the kinds of conditions end in the list of a system; the results of
// the elimination of the iso lines follow behind `isoLines`.
struct ConstraintSections {
    int anchors = 0;    // end of the anchors
    int interior = 0;   // end of the transitions inside the chains
    int chainEnds = 0;  // end of the eliminated transitions at the chain ends
    int isoLines = 0;   // end of the conditions "same coordinate as the reference of the curve"
};

// The system of the parameterisation of one part: two variables (u, v) per
// point of the cut copy and two offsets per chain along a disc cut; the
// conditions that tie them together; the energy as weighted rows over the
// variables that stay unknown; and which unknowns have to become whole
// numbers. The solving step reads it and changes weights and the normal
// equations.
struct System {
    // Variable 2w is u and 2w + 1 is v of copy point w; behind them follow
    // the two offsets of every chain along a disc cut.
    int copyPointCount = 0;
    int variableCount = 0;

    // Per variable its kind (variable_flag) and its grid period: a variable
    // with period m is rounded to multiples of m (0: continuous).
    std::vector<std::uint8_t> flags;
    std::vector<int> periods;
    // How many periods the conditions changed from 1 to something larger.
    int periodChanges = 0;

    ConstraintSet constraints;
    ConstraintSections sections;
    // Whether the conditions at the chain ends and those of the iso lines
    // could be eliminated (false: a pivot was too small; the conditions that
    // were not reached stand in the list without pivot and without terms).
    bool chainEndsEliminated = true;
    bool isoLinesEliminated = true;

    // Per curve chain the variable all its points are tied to, -1 when none
    // of its variables was free.
    std::vector<int> curveReference;

    // Energy: row i asks for sum of coefficient * unknown = targets[i] with
    // weights[i]; every row is sorted by unknown. Four rows per computed
    // triangle (triangleRowStart: triangle -> its first row, one entry more
    // than triangles), behind them the rows of iso lines that could not be
    // eliminated, with a very large weight.
    solver::ListMatrix rows;
    std::vector<double> targets;
    std::vector<double> weights;
    std::vector<int> triangleRowStart;
    int energyRowCount = 0;
    int penaltyRowCount = 0;

    // Numbering of the unknowns: variables that are neither eliminated nor
    // lie in the mirrored half, in ascending order (-1: no unknown).
    std::vector<int> variableToUnknown;
    std::vector<int> unknownToVariable;
    // Terms of rows whose variable has no unknown. A system that is handed
    // on has none.
    int unresolvedTerms = 0;

    // Normal equations of the rows, with the transposed rows they were built from.
    solver::NormalEquations normal;

    // Unknowns that have to become whole numbers, ascending, and their grid
    // periods (empty unless a period was changed).
    std::vector<int> integerUnknowns;
    std::vector<double> integerPeriods;

    [[nodiscard]] int unknownCount() const noexcept { return static_cast<int>(unknownToVariable.size()); }
};

} // namespace areko::remesher::param
