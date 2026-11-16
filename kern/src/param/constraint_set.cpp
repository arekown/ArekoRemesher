#include "constraint_set.h"

#include "param/constraint_terms.h"

#include <utility>

namespace areko::remesher::param {

TermList& ConstraintSet::add(int pivot, TermList terms)
{
    ofVariable_[static_cast<std::size_t>(pivot)] = count();
    list_.push_back({pivot, std::move(terms)});
    return list_.back().terms;
}

void ConstraintSet::append(std::vector<Constraint> conditions)
{
    for (Constraint& condition : conditions) {
        list_.push_back(std::move(condition));
    }
}

void ConstraintSet::removeAnchorTerms(TermList& terms) const
{
    for (std::size_t i = terms.size(); i-- > 0;) {
        const int condition = ofVariable(terms[i].index);
        if (condition >= 0 && condition < anchorCount_) {
            terms.erase(terms.begin() + static_cast<std::ptrdiff_t>(i));
        }
    }
}

void ConstraintSet::substituted(TermList& terms, std::span<const Constraint> sources, int offset, int first) const
{
    const int sourceEnd = static_cast<int>(sources.size());
    const int variableEnd = static_cast<int>(ofVariable_.size());
    // Only the terms the list starts with are candidates.
    std::size_t candidates = terms.size();
    std::size_t i = 0;
    while (i < candidates) {
        const Term term = terms[i];
        const int source = term.index < variableEnd ? ofVariable(term.index) - offset : -1;
        if (source < first || source >= sourceEnd) {
            ++i;
            continue;
        }
        terms.erase(terms.begin() + static_cast<std::ptrdiff_t>(i));
        --candidates;
        const TermList& replacement = sources[static_cast<std::size_t>(source)].terms;
        const std::size_t count = replacement.size();
        for (std::size_t k = 0; k < count; ++k) {
            const Term part = replacement[k];
            terms.push_back({part.index, part.value * term.value});
        }
    }
    ConstraintTerms::normalize(terms);
}

void ConstraintSet::resolve(int first, int end)
{
    for (int k = 0; k < end; ++k) {
        substituted(list_[static_cast<std::size_t>(k)].terms, list_, 0, first);
    }
}

void ConstraintSet::resolveAmong(int first)
{
    for (int k = count(); k-- > first;) {
        substituted(list_[static_cast<std::size_t>(k)].terms, list_, 0, first);
    }
}

void ConstraintSet::substituteInto(std::span<TermList> rows) const
{
    for (TermList& row : rows) {
        substituted(row, list_, 0, 0);
    }
}

} // namespace areko::remesher::param
