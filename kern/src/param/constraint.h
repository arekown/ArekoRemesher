#pragma once

#include "solver/sparse_matrix.h"

#include <vector>

namespace areko::remesher::param {

// One term of a linear expression: `index` is the variable, `value` its coefficient.
using Term = solver::SparseEntry;

// A linear expression. As a row it means "sum of value * x[index] = 0".
using TermList = solver::ListMatrix::List;

// An eliminated variable written as a sum over other variables:
// x[pivot] = sum of value * x[index], without a constant. pivot -1 = empty.
struct Constraint {
    int pivot = -1;
    TermList terms;
};

} // namespace areko::remesher::param
