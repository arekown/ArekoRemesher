#pragma once

#include "ngon_split_templates.h"

#include <vector>

namespace areko::remesher::extract {

// Rules for a polygon whose corners are four corners and straights only.
class StripGrid {
public:
    // The first two neighbouring corner pairs split the outline into two
    // sides; a chord between them closes a strip (or three quads for sides
    // of 3 and 5 edges). False when the rule does not apply.
    static bool strip(NgonSplitWork& work, int p);
    // Opposite sides of equal length: the polygon becomes a smoothed grid of
    // quads. False when the sides do not match.
    static bool grid(NgonSplitWork& work, int p);

private:
    static bool pairSides(NgonSplitWork& work, int p, int a, int b, int d1, int d2);
    // Arc length from local corner `from` to local corner `to` (round the end),
    // as share of the whole arc; `to` gets 1.
    static void arc(const NgonSplitWork& work, int p, int from, int to, std::vector<float>& share);
    static void smoothGrid(NgonSplitWork& work, const std::vector<int>& grid, int rows, int columns);
};

} // namespace areko::remesher::extract
