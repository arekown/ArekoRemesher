#pragma once

#include "field/cross_field.h"

#include <array>
#include <vector>

namespace areko::remesher::param {

// The parameterisation of one part: the solved system as grid coordinates on
// the cut copy. The later quads are the unit cells of this grid. The
// extraction reads it.
struct Solution {
    // (u, v) per point of the cut copy. Unknowns that had to become whole
    // numbers are whole. With a mirrored half its corners carry the
    // coordinates of their mirror corners with u negated.
    std::vector<std::array<float, 2>> uv;

    // Singular points of the part mesh the extraction starts from: those of
    // the chains; on a symmetric part the ones that touch the computed half,
    // each followed by its mirror point (then no longer ascending).
    field::Singularities singularities;
    // Further start points of the extraction: the points where a run of
    // curve edges changes its direction; without any singular point also the
    // border points with whole coordinates; on a symmetric part mirrored like
    // the singular points.
    std::vector<int> seeds;

    // Rounds of the driver, and the unknowns of its integer list: each was
    // fixed, on a whole number or, along a frozen border, on the value the
    // rounding of the borders kept for it.
    int rounds = 0;
    int integerUnknowns = 0;
};

} // namespace areko::remesher::param
