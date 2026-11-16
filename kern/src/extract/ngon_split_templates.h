#pragma once

#include "ngon_new_points.h"
#include "ngon_polygon_list.h"

#include <initializer_list>
#include <vector>

namespace areko::remesher::extract {

// The polygon list being split, its corner data and the new points.
struct NgonSplitWork {
    NgonPolygonList& list;
    NgonCornerData& corners;
    NgonNewPoints& points;
};

// Replacements of one polygon p of the list. Each removes p (the later
// polygons move forward) and appends its parts behind all others. Positions
// are absolute positions of the corner list unless said otherwise. No
// replacement adds corner data for its parts.
class SplitTemplates {
public:
    // Chord between the local corners a and b (0 .. n-1, any order).
    static void cut(NgonSplitWork& work, int p, int a, int b);
    // Corner b between a and c: a new point opposite b closes the quad
    // (a, b, c, X); the rest is split at the corner half way round.
    static void singleCorner(NgonSplitWork& work, int p, int a, int b, int c);
    // Three corners v1 v2 v3 between the straights v0 and v4: a new point
    // inside closes the quads (v0, v1, v2, X) and (v2, v3, v4, X).
    static void threeCorners(NgonSplitWork& work, int p, int v0, int v1, int v2, int v3, int v4);
    // Corners k and m with straights between them: two new points cut off
    // the quads at k and m.
    static void cornerStraightsCorner(NgonSplitWork& work, int p, int k, int m);
    // Three parts between the corners a, b and c (no new points).
    static void fan(NgonSplitWork& work, int p, int a, int b, int c);
    // Triangles of the planar tessellation, else a fan from the first corner.
    static void triangulate(NgonSplitWork& work, int p);
    // Removes polygon p and its corner data, then appends the parts.
    static void replace(NgonSplitWork& work, int p, const std::vector<std::vector<int>>& parts);
    // Point ids from position `from` round to position `to` (both included).
    static void walk(const NgonPolygonList& list, int p, int from, int to, std::vector<int>& out);
};

} // namespace areko::remesher::extract
