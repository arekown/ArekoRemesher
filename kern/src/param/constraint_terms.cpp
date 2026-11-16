#include "constraint_terms.h"

#include "util/rounding.h"

#include <cmath>
#include <utility>

namespace areko::remesher::param {

void ConstraintTerms::normalize(TermList& terms)
{
    solver::ListMatrix::order(terms, solver::ListMatrix::SameIndex::SummedNonZero);
}

void ConstraintTerms::dropTiny(TermList& terms, double limit)
{
    bool found = false;
    for (Term& term : terms) {
        if (limit > std::fabs(term.value) && term.value != 0.0) {
            term.value = 0.0;
            found = true;
        }
    }
    if (found) {
        std::erase_if(terms, [](const Term& term) { return term.value == 0.0; });
    }
}

TermList ConstraintTerms::combine(double factorFirst, const TermList& first, double factorSecond,
                                  const TermList& second)
{
    TermList out;
    out.reserve(first.size() + second.size());
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < first.size() && j < second.size()) {
        const Term& a = first[i];
        const Term& b = second[j];
        if (a.index == b.index) {
            const double sum = a.value * factorFirst + b.value * factorSecond;
            if (sum != 0.0) {
                out.push_back({a.index, sum});
            }
            ++i;
            ++j;
        } else if (a.index < b.index) {
            out.push_back({a.index, a.value * factorFirst});
            ++i;
        } else {
            out.push_back({b.index, b.value * factorSecond});
            ++j;
        }
    }
    for (; i < first.size(); ++i) {
        out.push_back({first[i].index, first[i].value * factorFirst});
    }
    for (; j < second.size(); ++j) {
        out.push_back({second[j].index, second[j].value * factorSecond});
    }
    return out;
}

void ConstraintTerms::reduceToIntegers(TermList& row)
{
    if (row.empty()) {
        return;
    }
    // A NaN magnitude never becomes the smallest one.
    double smallest = magnitudeStart;
    for (const Term& term : row) {
        const double magnitude = std::fabs(term.value);
        if (magnitude == 1.0) {
            return;
        }
        if (magnitude < smallest) {
            smallest = magnitude;
        }
    }
    TermList reduced = row;
    for (Term& term : reduced) {
        const double quotient = term.value / smallest;
        if (static_cast<double>(util::truncateToInt(quotient)) != quotient) {
            return;
        }
        term.value = quotient;
    }
    row = std::move(reduced);
}

} // namespace areko::remesher::param
