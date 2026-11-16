#include "transition_constraints.h"

#include "field/curve_alignment.h"
#include "mesh/vec3.h"
#include "param/constraint_terms.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <utility>

namespace areko::remesher::param {

TransitionConstraints::Turn TransitionConstraints::turnOf(float rotation) noexcept
{
    // Truncated the way the processor does it: a value that is no number or
    // out of range gives the smallest int, which counts as no turn. The
    // remainder is taken without sign, so that nothing overflows.
    const std::int32_t r = mesh::truncToInt(rotation);
    const auto bits = static_cast<std::uint32_t>(r);
    const std::uint32_t quarter = (r >= 0 ? bits : 0u - 3u * bits) & 3u;
    static constexpr std::array<Turn, 4> turns{
        {{false, 1.0, 0.0}, {true, 0.0, 1.0}, {false, -1.0, 0.0}, {true, 0.0, -1.0}}};
    return turns[quarter];
}

TransitionConstraints::Sides TransitionConstraints::sidesOf(int chain) const
{
    const Chains& chains = context_.chains();
    const std::vector<int>& a = chains.sideA[static_cast<std::size_t>(chain)];
    const std::vector<int>& b = chains.sideB[static_cast<std::size_t>(chain)];
    if (a.empty() || b.empty() || a.size() != b.size()) {
        return {};
    }
    return {&a, &b};
}

bool TransitionConstraints::anyMirrored(std::initializer_list<int> copyPoints) const
{
    return context_.symmetric() &&
           std::any_of(copyPoints.begin(), copyPoints.end(), [&](int w) { return context_.mirrored(w); });
}

Term TransitionConstraints::turned(const Turn& turn, int coordinate, int b, double sign) noexcept
{
    if (coordinate == 0) {
        return turn.odd ? Term{2 * b + 1, sign * -turn.sine} : Term{2 * b, sign * turn.cosine};
    }
    return turn.odd ? Term{2 * b, sign * turn.sine} : Term{2 * b + 1, sign * turn.cosine};
}

// The one place that orders the terms of a transition, whatever the turn:
// equal variables are added up and terms that cancel are gone. A row may end
// up without terms; it stays a row all the same.
void TransitionConstraints::finish(TermList& row) const
{
    ConstraintTerms::normalize(row);
    system_.constraints.removeAnchorTerms(row);
}

void TransitionConstraints::addInterior()
{
    const Chains& chains = context_.chains();
    for (int chain = 0; chain < chains.count(); ++chain) {
        const std::size_t c = static_cast<std::size_t>(chain);
        if (chains.fullTurnChains.testBounded(c) || chains.curveCheckChains.testBounded(c)) {
            continue;
        }
        const Sides sides = sidesOf(chain);
        if (!sides.usable()) {
            continue;
        }
        const std::vector<int>& a = *sides.a;
        const std::vector<int>& b = *sides.b;
        if (anyMirrored({a.front(), b.front()})) {
            continue;
        }
        const Turn turn = turnOf(chains.rotation[c]);
        for (std::size_t k = 1; k + 1 < a.size(); ++k) {
            if (anyMirrored({a[k], b[k]})) {
                continue;
            }
            for (int coordinate = 0; coordinate < 2; ++coordinate) {
                // The condition is entered for its pivot before the terms of
                // anchors are dropped: a pivot that was an anchor is none
                // any more, and a term that names it stays.
                TermList& terms = system_.constraints.add(2 * a[k] + coordinate);
                terms = {turned(turn, coordinate, b[k], 1.0), turned(turn, coordinate, b.front(), -1.0),
                         Term{2 * a.front() + coordinate, 1.0}};
                finish(terms);
            }
        }
    }
    system_.sections.interior = system_.constraints.count();
}

std::vector<TermList> TransitionConstraints::endRows() const
{
    const Chains& chains = context_.chains();
    const mesh::Adjacency& adjacency = context_.adjacency();
    std::vector<TermList> rows;
    // Every row is kept, also one whose terms cancel or are all anchors: the
    // elimination counts its rows, and an empty one gets no pivot.
    auto append = [&](TermList row) {
        finish(row);
        rows.push_back(std::move(row));
    };
    for (int chain = 0; chain < chains.count(); ++chain) {
        const std::size_t c = static_cast<std::size_t>(chain);
        if (chains.fullTurnChains.testBounded(c)) {
            continue;
        }
        const Sides sides = sidesOf(chain);
        if (!sides.usable()) {
            continue;
        }
        const int a0 = sides.a->front();
        const int a1 = sides.a->back();
        const int b0 = sides.b->front();
        const int b1 = sides.b->back();
        if (anyMirrored({a0, a1, b0, b1})) {
            continue;
        }
        const Turn turn = turnOf(chains.rotation[c]);
        const int offset = chains.offsetVariable[c];
        if (offset == 0) {
            // A chain of a curve check cut keeps only the coordinate that
            // does not change along its edge.
            bool rowFor[2] = {true, true};
            if (chains.curveCheckChains.testBounded(c)) {
                const int corner = adjacency.edgeCorners(chains.edges.row(chain)[0])[0];
                const bool alongFirst = field::CurveAlignment::followsFirst(context_.mesh(), *context_.part().crossField,
                                                                            context_.graph(), corner);
                rowFor[alongFirst ? 1 : 0] = false;
            }
            for (int coordinate = 0; coordinate < 2; ++coordinate) {
                if (rowFor[coordinate]) {
                    append({turned(turn, coordinate, b1, 1.0), turned(turn, coordinate, b0, -1.0),
                            Term{2 * a0 + coordinate, 1.0}, Term{2 * a1 + coordinate, -1.0}});
                }
            }
        } else {
            for (const auto& [a, b] : {std::pair{a0, b0}, std::pair{a1, b1}}) {
                for (int coordinate = 0; coordinate < 2; ++coordinate) {
                    append({turned(turn, coordinate, b, 1.0), Term{2 * a + coordinate, -1.0},
                            Term{offset + coordinate, 1.0}});
                }
            }
        }
    }
    return rows;
}

} // namespace areko::remesher::param
