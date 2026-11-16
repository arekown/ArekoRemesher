#pragma once

#include "ngon_corner_classifier.h"
#include "ngon_split_templates.h"

namespace areko::remesher::extract {

// Splits the polygons of a list until none has five or more corners. Each
// round applies the first rule that fits the corner flags of polygon p; its
// parts go behind all others and are classified, and p then names the next
// polygon.
class PolygonSplitter {
public:
    static void split(NgonPolygonList& list, NgonNewPoints& points);

private:
    static void choose(NgonSplitWork& work, int p, const FlagCounts& k);
    // Rules in their order of trial (true = applied).
    static bool cutOffQuad(NgonSplitWork& work, int p);
    static bool threeBetweenStraights(NgonSplitWork& work, int p);
    static bool cornerStraightsCorner(NgonSplitWork& work, int p);
    static bool singleCorner(NgonSplitWork& work, int p);
    // Cheapest chord i -> i+3, else the planar triangulation.
    static void bestChord(NgonSplitWork& work, int p);
    // Cheapest chord from a concave corner, else the planar triangulation.
    static void concaveChord(NgonSplitWork& work, int p);
};

} // namespace areko::remesher::extract
