#include "border_rounding.h"

#include "mesh/boundary_loops.h"
#include "util/rounding.h"

#include <cassert>
#include <cmath>
#include <initializer_list>
#include <span>
#include <vector>

namespace areko::remesher::param {

using namespace areko::remesher::mesh;

namespace {

// One coordinate of one copy of a kept border point.
struct Entry {
    double whole = 0.0; // the whole number nearest to its value
    double rest = 0.0;  // distance of its value from that number
    int unknown = -1;   // -1: a condition gives the coordinate
    int copyPoint = 0;
};

// Whether an edge of the cut copy joins two copy points. A point counts as
// joined with itself as soon as it has an edge.
[[nodiscard]] bool joined(const Adjacency& copy, int a, int b)
{
    return a == b ? !copy.pointEdges(a).empty() : copy.edgeBetween(a, b) != -1;
}

// The border rounding of one round.
class BorderSelection {
public:
    BorderSelection(const SolveState& state, const std::vector<FreezeMark>& marks, solver::DriverRequest& request)
        : state_(state),
          marks_(marks),
          request_(request),
          x_(*request.x),
          position_(static_cast<std::size_t>(state.system.unknownCount()), -1),
          chosen_(state.driver.integerUnknowns().size())
    {
        const std::vector<int>& unknowns = state.driver.integerUnknowns();
        for (std::size_t q = 0; q < unknowns.size(); ++q) {
            position_[static_cast<std::size_t>(unknowns[q])] = static_cast<int>(q);
        }
    }

    void run()
    {
        const BoundaryLoops loops(state_.context.mesh());
        for (int loop = 0; loop < loops.count(); ++loop) {
            for (const int coordinate : {0, 1}) {
                const std::vector<Entry> list = entries(loops.points().row(loop), coordinate);
                if (!list.empty()) {
                    select(list, keptAsIs(list, coordinate));
                }
            }
        }
    }

private:
    // The coordinate of every copy of the kept points of a loop, in walking order.
    [[nodiscard]] std::vector<Entry> entries(std::span<const int> loop, int coordinate) const
    {
        const System& system = state_.system;
        std::vector<Entry> list;
        for (const int p : loop) {
            if (marks_[static_cast<std::size_t>(p)].value != 1) {
                continue;
            }
            for (const int w : state_.context.copies(p)) {
                const int variable = 2 * w + coordinate;
                Entry entry;
                entry.copyPoint = w;
                entry.unknown = system.variableToUnknown[static_cast<std::size_t>(variable)];
                double value = 0.0;
                if (entry.unknown != -1) {
                    value = x_[static_cast<std::size_t>(entry.unknown)];
                } else {
                    const int condition = system.constraints.ofVariable(variable);
                    if (condition == -1) {
                        continue;
                    }
                    for (const Term& term : system.constraints.list()[static_cast<std::size_t>(condition)].terms) {
                        // A term over a variable without unknown adds nothing.
                        const int unknown = system.variableToUnknown[static_cast<std::size_t>(term.index)];
                        if (unknown != -1) {
                            value = term.value * x_[static_cast<std::size_t>(unknown)] + value;
                        }
                    }
                }
                // The whole number goes straight to double.
                entry.whole = static_cast<double>(util::roundHalfAway(value));
                entry.rest = std::fabs(value - entry.whole);
                list.push_back(entry);
            }
        }
        return list;
    }

    // The entries that keep their value as it is: one of every two neighbours
    // in conflict. Passes are repeated until one marks nothing new. Which
    // entry a conflict picks depends on the whole numbers, kinds, unknowns
    // and rests, never on the marks: the second pass finds nothing new, so
    // the loop ends after two passes at the most.
    [[nodiscard]] BitArray keptAsIs(const std::vector<Entry>& list, int coordinate) const
    {
        const Adjacency& copy = state_.context.copy().adjacency();
        auto frozenUnknown = [&](const Entry& entry) { return frozen(entry, coordinate) && entry.unknown != -1; };
        const std::size_t n = list.size();
        BitArray kept(n);
        bool marked = true;
        while (marked) {
            marked = false;
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t j = (i + 1) % n;
                const Entry& first = list[i];
                const Entry& second = list[j];
                if (first.whole != second.whole || !joined(copy, first.copyPoint, second.copyPoint)) {
                    continue;
                }
                std::size_t pick = j;
                if (frozenUnknown(first)) {
                    pick = frozenUnknown(second) && second.rest > first.rest ? j : i;
                } else if (!frozen(second, coordinate)) {
                    continue;
                }
                if (!kept.test(pick)) {
                    kept.set(pick);
                    marked = true;
                }
            }
        }
        return kept;
    }

    void select(const std::vector<Entry>& list, const BitArray& kept)
    {
        const BitArray& fixed = state_.driver.fixedPositions();
        for (std::size_t e = 0; e < list.size(); ++e) {
            const Entry& entry = list[e];
            if (entry.unknown == -1) {
                continue;
            }
            const int position = position_[static_cast<std::size_t>(entry.unknown)];
            // An unknown outside the integer list has no position to select.
            if (position == -1 || fixed.testBounded(static_cast<std::size_t>(position))) {
                continue;
            }
            // A point in two border loops is selected in the first one only.
            if (chosen_.test(static_cast<std::size_t>(position))) {
                continue;
            }
            chosen_.set(static_cast<std::size_t>(position));
            const double value = kept.test(e) ? x_[static_cast<std::size_t>(entry.unknown)] : entry.whole;
            request_.selection.push_back({position, value});
        }
    }

    [[nodiscard]] bool frozen(const Entry& entry, int coordinate) const
    {
        return state_.system.flags[static_cast<std::size_t>(2 * entry.copyPoint + coordinate)] ==
               variable_flag::frozenBorder;
    }

    const SolveState& state_;
    const std::vector<FreezeMark>& marks_;
    solver::DriverRequest& request_;
    const std::vector<double>& x_;
    // Unknown -> its position in the integer list, -1: none.
    std::vector<int> position_;
    // Positions selected in this round.
    BitArray chosen_;
};

} // namespace

RoundingStep::Outcome BorderRounding::select(const SolveState& state, solver::DriverRequest& request)
{
    // Only a part with frozen borders is asked.
    const std::vector<FreezeMark>* marks = state.context.part().frozenBorderMarks();
    assert(marks != nullptr);
    const std::size_t before = request.selection.size();
    BorderSelection(state, *marks, request).run();
    return request.selection.size() > before ? Outcome::Delivered : Outcome::Failed;
}

} // namespace areko::remesher::param
