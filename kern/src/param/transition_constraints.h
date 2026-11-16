#pragma once

#include "param/system.h"
#include "param/system_context.h"

#include <initializer_list>
#include <vector>

namespace areko::remesher::param {

// Transitions across the cut chains: the coordinates on side A of a chain
// are those on side B turned by the quarter turns of the chain,
//     uv_A(k) - uv_A(0) = R * (uv_B(k) - uv_B(0)),   R = [[cos, -sin], [sin, cos]],
// with a free shift between the sides. A chain along a disc cut carries the
// shift as its two offset variables instead: uv_A = R * uv_B + offset.
class TransitionConstraints {
public:
    TransitionConstraints(const SystemContext& context, System& system) : context_(context), system_(system) {}

    // Cosine and sine of a number of quarter turns: whole numbers.
    struct Turn {
        bool odd = false;
        double cosine = 1.0;
        double sine = 0.0;
    };
    // The whole part of `rotation` counts; a value that is no number or too
    // large for a whole number counts as no turn.
    [[nodiscard]] static Turn turnOf(float rotation) noexcept;

    // Inside the chains: for every chain point between the ends one
    // condition per coordinate with the point on side A as pivot. Chains at
    // a whole turn and chains of a curve check cut get none, nor do chains
    // without matching sides or with a mirrored first point. Conditions are
    // not put into each other; only terms of anchors are dropped, after the
    // condition is entered for its pivot.
    void addInterior();

    // At the chain ends: rows "= 0" for the elimination, in chain order.
    //  - a chain without offset ties its last point to its first (one row per
    //    coordinate; a chain of a curve check cut only the coordinate that
    //    stays the same across its edge),
    //  - a chain with offsets ties each of its two ends across the cut
    //    (two rows per end).
    // Chains at a whole turn get none, nor do chains without matching sides
    // or with a mirrored end point. A row without terms stays in the list.
    [[nodiscard]] std::vector<TermList> endRows() const;

private:
    // The two sides of a chain when it takes part: both present and equally long.
    struct Sides {
        const std::vector<int>* a = nullptr;
        const std::vector<int>* b = nullptr;
        [[nodiscard]] bool usable() const noexcept { return a != nullptr; }
    };
    [[nodiscard]] Sides sidesOf(int chain) const;
    [[nodiscard]] bool anyMirrored(std::initializer_list<int> copyPoints) const;
    // Coordinate `coordinate` (0 or 1) of R * uv(b) as one term: the turn
    // decides which coordinate of b it reads and with which sign.
    [[nodiscard]] static Term turned(const Turn& turn, int coordinate, int b, double sign) noexcept;
    // Brings a row into its normal form and drops the terms of anchors.
    void finish(TermList& row) const;

    const SystemContext& context_;
    System& system_;
};

} // namespace areko::remesher::param
