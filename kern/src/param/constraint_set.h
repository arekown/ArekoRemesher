#pragma once

#include "param/constraint.h"

#include <span>
#include <vector>

namespace areko::remesher::param {

// The conditions of a system: eliminated variables written as sums over other
// variables, in the order they were found, with the map from every variable
// to the condition it is the pivot of (-1: none; the last entry made counts).
// The first conditions are the anchors: variables set to zero.
class ConstraintSet {
public:
    ConstraintSet() = default;
    explicit ConstraintSet(int variableCount) : ofVariable_(static_cast<std::size_t>(variableCount), -1) {}

    [[nodiscard]] int count() const noexcept { return static_cast<int>(list_.size()); }
    [[nodiscard]] const std::vector<Constraint>& list() const noexcept { return list_; }
    [[nodiscard]] int ofVariable(int variable) const { return ofVariable_[static_cast<std::size_t>(variable)]; }
    // The map itself, for an elimination that enters its pivots.
    [[nodiscard]] std::vector<int>& variableMap() noexcept { return ofVariable_; }
    [[nodiscard]] int anchorCount() const noexcept { return anchorCount_; }

    // Appends a condition and makes it the condition of its pivot. Returns
    // its terms, for a caller that fills them in after the entry; they stay
    // in place until the next condition is appended.
    TermList& add(int pivot, TermList terms = {});
    // Appends the conditions of an elimination, whose pivots it has entered
    // in the map.
    void append(std::vector<Constraint> conditions);
    // All conditions so far are anchors.
    void closeAnchors() noexcept { anchorCount_ = count(); }

    // Removes the terms of anchor variables; the order of the rest stays.
    void removeAnchorTerms(TermList& terms) const;

    // Puts conditions into a term list: every term of the list as it stands
    // whose variable is the pivot of a condition in `sources` is replaced by
    // the terms of that condition times its coefficient, appended at the end;
    // appended terms are not looked at again. The list is normalised after
    // that. The condition of a variable is its map entry minus `offset`,
    // taken as an index into `sources`; it counts from `first` on. A term
    // whose variable lies beyond the map stays as it is, in rows (which may
    // carry a column of constants there) and in conditions alike.
    void substituted(TermList& terms, std::span<const Constraint> sources, int offset, int first) const;

    // Puts the conditions from `first` on into the conditions [0, end).
    void resolve(int first, int end);
    // Puts the conditions from `first` on into each other, from the last
    // one backwards.
    void resolveAmong(int first);
    // Puts all conditions into every row.
    void substituteInto(std::span<TermList> rows) const;

private:
    std::vector<Constraint> list_;
    std::vector<int> ofVariable_;
    int anchorCount_ = 0;
};

} // namespace areko::remesher::param
