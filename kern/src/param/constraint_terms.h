#pragma once

#include "param/constraint.h"

#include <limits>

namespace areko::remesher::param {

// Arithmetic on term lists. Comparisons with 0 and 1 are exact; a NaN
// coefficient is neither 0 nor 1 and is kept wherever zeros are dropped.
class ConstraintTerms {
public:
    // Start of every search for the smallest magnitude of a coefficient.
    static constexpr double magnitudeStart = static_cast<double>(std::numeric_limits<float>::max());

    // Brings a term list into its normal form: sorted by variable, terms of
    // the same variable added up in sorted order, coefficients equal to 0
    // skipped and sums equal to 0 removed. This is the one way to order the
    // terms of a condition, whatever they were built from.
    static void normalize(TermList& terms);

    // Removes every term whose coefficient is not 0 but smaller in magnitude
    // than `limit`; if there was one, all terms that are exactly 0 go with
    // them. The order of the rest is kept.
    static void dropTiny(TermList& terms, double limit);

    // first * factorFirst + second * factorSecond, both sorted by variable. A
    // variable of both lists is dropped when its sum is exactly 0; a variable
    // of one list only is kept even with a product of 0.
    [[nodiscard]] static TermList combine(double factorFirst, const TermList& first, double factorSecond,
                                          const TermList& second);

    // Divides all coefficients by the smallest magnitude if every quotient is
    // a whole number that fits an int. A row that holds a coefficient of
    // magnitude 1 stays as it is.
    static void reduceToIntegers(TermList& row);
};

} // namespace areko::remesher::param
