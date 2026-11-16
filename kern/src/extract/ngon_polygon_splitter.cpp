#include "ngon_polygon_splitter.h"

#include "ngon_cut_check.h"
#include "ngon_strip_grid.h"

#include <limits>

namespace areko::remesher::extract {

namespace {

constexpr float cornerPenalty = 1.5f;
constexpr float concaveBonus = 0.5f;

} // namespace

void PolygonSplitter::split(NgonPolygonList& list, NgonNewPoints& points)
{
    NgonCornerData corners;
    CornerClassifier::classify(list, points, corners);
    NgonSplitWork work{list, corners, points};
    int polygons = list.polygonCount();
    for (int p = 0; p < polygons; ++p) {
        while (list.size(p) >= 5) {
            FlagCounts k = CornerClassifier::count(corners, list.begin(p), list.end(p));
            if (k.unclear != 0) {
                CornerClassifier::reclassify(corners, k, list.begin(p), list.end(p));
            }
            choose(work, p, k);
            CornerClassifier::classify(list, points, corners);
            polygons = list.polygonCount();
        }
    }
}

void PolygonSplitter::choose(NgonSplitWork& work, int p, const FlagCounts& k)
{
    const NgonPolygonList& list = work.list;
    const int s = list.begin(p);
    const int n = list.size(p);
    if (k.concave != 0) {
        concaveChord(work, p);
        return;
    }
    if (n == 5) {
        if (k.corners == 4 && k.straights == 1) {
            SplitTemplates::cut(work, p, k.lastStraight - s, list.wrap(p, k.lastStraight + 3) - s);
            return;
        }
        bestChord(work, p);
        return;
    }
    if (k.corners == 4 && n == k.straights + 4) {
        if (StripGrid::strip(work, p) || StripGrid::grid(work, p)) {
            return;
        }
    }
    if (cutOffQuad(work, p) || threeBetweenStraights(work, p) || cornerStraightsCorner(work, p) ||
        singleCorner(work, p)) {
        return;
    }
    bestChord(work, p);
}

bool PolygonSplitter::cutOffQuad(NgonSplitWork& work, int p)
{
    const NgonPolygonList& list = work.list;
    const auto& flag = work.corners.flag;
    auto f = [&](int i) { return flag[static_cast<std::size_t>(i)]; };
    auto cornerLike = [&](int i) { return f(i) == 1 || f(i) == 2; };
    const int s = list.begin(p);
    for (int k = s; k < list.end(p); ++k) {
        if (cornerLike(k)) {
            continue;
        }
        const int a = list.next(p, k);
        const int b = list.next(p, a);
        const int c = list.next(p, b);
        if (cornerLike(c)) {
            continue;
        }
        if (f(a) == 1 && f(b) == 1 && (f(k) != 0 || f(c) != 0)) {
            SplitTemplates::cut(work, p, k - s, c - s);
            return true;
        }
    }
    return false;
}

bool PolygonSplitter::threeBetweenStraights(NgonSplitWork& work, int p)
{
    const NgonPolygonList& list = work.list;
    const auto& flag = work.corners.flag;
    auto f = [&](int i) { return flag[static_cast<std::size_t>(i)]; };
    for (int k = list.begin(p); k < list.end(p); ++k) {
        if (f(k) != 0) {
            continue;
        }
        const int a = list.next(p, k);
        const int b = list.next(p, a);
        const int c = list.next(p, b);
        const int d = list.next(p, c);
        if (f(a) == 1 && f(b) == 1 && f(c) == 1 && f(d) == 0) {
            SplitTemplates::threeCorners(work, p, k, a, b, c, d);
            return true;
        }
    }
    return false;
}

bool PolygonSplitter::cornerStraightsCorner(NgonSplitWork& work, int p)
{
    const NgonPolygonList& list = work.list;
    const auto& flag = work.corners.flag;
    auto f = [&](int i) { return flag[static_cast<std::size_t>(i)]; };
    const int n = list.size(p);
    for (int k = list.begin(p); k < list.end(p); ++k) {
        if (f(k) != 1) {
            continue;
        }
        const int stop = list.wrap(p, k + n - 2);
        int m = list.next(p, k);
        if (m == stop || f(m) != 0) {
            continue; // a corner right behind k leaves no straights between
        }
        int straights = 0;
        while (true) {
            ++straights;
            m = list.next(p, m);
            if (m == stop || f(m) != 0) {
                break;
            }
        }
        if (m != stop && f(m) == 1 && straights >= 2) {
            SplitTemplates::cornerStraightsCorner(work, p, k, m);
            return true;
        }
    }
    return false;
}

bool PolygonSplitter::singleCorner(NgonSplitWork& work, int p)
{
    const NgonPolygonList& list = work.list;
    const auto& flag = work.corners.flag;
    auto f = [&](int i) { return flag[static_cast<std::size_t>(i)]; };
    for (int k = list.begin(p); k < list.end(p); ++k) {
        if (f(k) != 1) {
            continue;
        }
        const int before = list.prev(p, k);
        const int after = list.next(p, k);
        if (f(before) == 0 && f(after) == 0) {
            SplitTemplates::singleCorner(work, p, before, k, after);
            return true;
        }
    }
    return false;
}

void PolygonSplitter::bestChord(NgonSplitWork& work, int p)
{
    const NgonPolygonList& list = work.list;
    const auto& flag = work.corners.flag;
    const int s = list.begin(p);
    float best = std::numeric_limits<float>::max();
    int bestAt = -1;
    for (int i = s; i < list.end(p); ++i) {
        const int j = list.wrap(p, i + 3);
        const CutCheck::Chord chord = CutCheck::chord(list, work.points, p, i, j);
        if (chord.blocked) {
            continue;
        }
        float cost = chord.length;
        if (flag[static_cast<std::size_t>(i)] == 1) {
            cost = cost * cornerPenalty;
        }
        if (flag[static_cast<std::size_t>(j)] == 1) {
            cost = cost * cornerPenalty;
        }
        cost = CutCheck::angleFactor(list, work.points, p, i, j, chord.length, cost);
        if (cost < best) {
            best = cost;
            bestAt = i;
        }
    }
    if (bestAt == -1) {
        SplitTemplates::triangulate(work, p);
        return;
    }
    SplitTemplates::cut(work, p, bestAt - s, list.wrap(p, bestAt + 3) - s);
}

void PolygonSplitter::concaveChord(NgonSplitWork& work, int p)
{
    const NgonPolygonList& list = work.list;
    const auto& flag = work.corners.flag;
    const int s = list.begin(p);
    const int n = list.size(p);
    float best = std::numeric_limits<float>::max();
    int bestI = -1;
    int bestJ = -1;
    for (int i = s; i < list.end(p); ++i) {
        if (flag[static_cast<std::size_t>(i)] != -1) {
            continue;
        }
        const int stop = list.wrap(p, i + n - 2);
        for (int j = list.wrap(p, i + 3); j != stop; j = list.next(p, j)) {
            const CutCheck::Chord chord = CutCheck::chord(list, work.points, p, i, j);
            if (chord.blocked) {
                continue;
            }
            float cost = chord.length;
            if (flag[static_cast<std::size_t>(j)] != -1) {
                cost = cost * concaveBonus;
            }
            if (cost < best) {
                best = cost;
                bestI = i;
                bestJ = j;
            }
        }
    }
    if (bestI == -1) {
        SplitTemplates::triangulate(work, p);
        return;
    }
    SplitTemplates::cut(work, p, bestI - s, bestJ - s);
}

} // namespace areko::remesher::extract
