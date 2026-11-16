#include "ngon_strip_grid.h"

#include <cmath>
#include <limits>
#include <utility>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

constexpr float evenPenaltyFactor = 1.4f;
constexpr float evenPenaltyOffset = 0.15f;
constexpr float mirrorBonus = 0.85f;
constexpr int gridIterations = 10;
// Weights of the grid smoothing along rows and columns; fixed to the middle,
// so edges of length zero do not matter.
constexpr float rowWeight = 0.5f;
constexpr float columnWeight = 0.5f;

} // namespace

bool StripGrid::strip(NgonSplitWork& work, int p)
{
    const NgonPolygonList& list = work.list;
    const auto& flag = work.corners.flag;
    const int s = list.begin(p);
    const int e = list.end(p);
    const int n = e - s;
    int a = -1;
    int k = s;
    while (k < e) {
        if (flag[static_cast<std::size_t>(k)] == 1 && flag[static_cast<std::size_t>(list.next(p, k))] == 1) {
            if (a == -1) {
                a = k;
                k = a + 2;
                continue;
            }
            const int b = k;
            const int d1 = b - a;
            const int d2 = n - d1;
            if (d1 <= 2 || d2 <= 2) {
                return false;
            }
            int x = -1;
            if (d2 == 5 && d1 == 3) {
                x = list.wrap(p, a + 2);
            } else if (d2 == 3 && d1 == 5) {
                x = list.wrap(p, b + 2);
            }
            if (x >= 0) {
                const int y = list.wrap(p, x + 3);
                SplitTemplates::fan(work, p, x, y, list.wrap(p, y + 2));
                return true;
            }
            return pairSides(work, p, a, b, d1, d2);
        }
        ++k;
    }
    return false;
}

void StripGrid::arc(const NgonSplitWork& work, int p, int from, int to, std::vector<float>& share)
{
    const NgonPolygonList& list = work.list;
    const int s = list.begin(p);
    const int n = list.size(p);
    if (from != to) {
        float sum = 0.0f;
        int i = from;
        int q = from;
        do {
            q = (i + 1) % n;
            share[static_cast<std::size_t>(i)] = sum;
            sum += distance(work.points.position(list.vertex(s + i)), work.points.position(list.vertex(s + q)));
            i = q;
        } while (q != to);
        for (i = from; i != to; i = (i + 1) % n) {
            share[static_cast<std::size_t>(i)] /= sum;
        }
    }
    share[static_cast<std::size_t>(to)] = 1.0f;
}

bool StripGrid::pairSides(NgonSplitWork& work, int p, int a, int b, int d1, int d2)
{
    const NgonPolygonList& list = work.list;
    const int s = list.begin(p);
    const int n = list.size(p);
    std::vector<float> share(static_cast<std::size_t>(n), 0.0f);
    arc(work, p, (a - s + 1) % n, b - s, share);
    arc(work, p, (b - s + 1) % n, a - s, share);
    auto at = [&](int position) { return share[static_cast<std::size_t>(position - s)]; };
    const bool mirrored = d1 == d2;
    float best = std::numeric_limits<float>::max();
    float worst = 0.0f;
    int bestI = -1;
    int bestJ = -1;
    int step = 1;
    for (int i = list.wrap(p, a + 2); i != b; i = list.next(p, i), ++step) {
        int partnerCountdown = d1 - 2;
        for (int j = list.wrap(p, b + 2); j != a; j = list.next(p, j), --partnerCountdown) {
            float c = std::fabs((1.0f - at(j)) - at(i));
            if (((((i - j) + n) | n) & 1) == 0) {
                c = c * evenPenaltyFactor + evenPenaltyOffset;
            }
            if (mirrored && partnerCountdown == step) {
                c = c * mirrorBonus;
                if (c > worst) {
                    worst = c;
                }
            }
            if (c < best) {
                best = c;
                bestI = i;
                bestJ = j;
            }
        }
    }
    if (bestI == -1 || (mirrored && 0.5f / static_cast<float>(d1 - 1) > worst)) {
        return false;
    }
    SplitTemplates::cut(work, p, bestI - s, bestJ - s);
    return true;
}

bool StripGrid::grid(NgonSplitWork& work, int p)
{
    const NgonPolygonList& list = work.list;
    const auto& flag = work.corners.flag;
    const int s = list.begin(p);
    const int e = list.end(p);
    const int n = e - s;
    std::vector<int> corners;
    int sides[4] = {0, 0, 0, 0};
    int count = 0;
    for (int q = s + 1;; ++q) {
        ++sides[count & 3];
        const int at = q >= e ? s : q;
        if (flag[static_cast<std::size_t>(at)] == 1) {
            ++count;
            corners.push_back(at);
        }
        if (q == e) {
            break;
        }
    }
    if (count != 4 || sides[0] != sides[2] || sides[1] != sides[3]) {
        return false;
    }
    const int across = sides[0]; // edges E4 -> E1
    const int down = sides[1];   // edges E1 -> E2
    const int columns = across + 1;
    const int rows = down + 1;
    std::vector<int> outline;
    outline.reserve(static_cast<std::size_t>(n));
    for (int i = corners[3]; i < e; ++i) {
        outline.push_back(list.vertex(i));
    }
    for (int i = s; i < corners[3]; ++i) {
        outline.push_back(list.vertex(i));
    }
    std::vector<int> g(static_cast<std::size_t>(rows * columns), -1);
    auto cell = [&](int r, int c) -> int& { return g[static_cast<std::size_t>(r * columns + c)]; };
    auto border = [&](int i) { return outline[static_cast<std::size_t>(i)]; };
    for (int c = 0; c <= across; ++c) {
        cell(0, c) = border(c);
        cell(down, c) = border(2 * across + down - c);
    }
    for (int r = 1; r < down; ++r) {
        cell(r, across) = border(across + r);
        cell(r, 0) = border(n - r);
    }
    NgonNewPoints& points = work.points;
    for (int r = 1; r < down; ++r) {
        const int left = cell(r, 0);
        const int right = cell(r, across);
        const Vec3f pl = points.position(left);
        const Vec3f pr = points.position(right);
        const float v = (static_cast<float>(r) * 0.5f) / static_cast<float>(down);
        const float vRest = 0.5f - v;
        for (int c = 1; c < across; ++c) {
            const int top = cell(0, c);
            const int bottom = cell(down, c);
            const float u = (static_cast<float>(c) * 0.5f) / static_cast<float>(across);
            const float uRest = 0.5f - u;
            const Vec3f x = ((pr * u + pl * uRest) + points.position(top) * vRest) + points.position(bottom) * v;
            cell(r, c) = points.add(x);
            points.addParent(left, uRest);
            points.addParent(right, u);
            points.addParent(top, vRest);
            points.addParent(bottom, v);
            points.endParents();
        }
    }
    smoothGrid(work, g, rows, columns);
    std::vector<std::vector<int>> quads;
    quads.reserve(static_cast<std::size_t>(down * across));
    for (int r = 0; r < down; ++r) {
        for (int c = 0; c < across; ++c) {
            quads.push_back({cell(r, c), cell(r, c + 1), cell(r + 1, c + 1), cell(r + 1, c)});
        }
    }
    SplitTemplates::replace(work, p, quads);
    return true;
}

void StripGrid::smoothGrid(NgonSplitWork& work, const std::vector<int>& grid, int rows, int columns)
{
    std::vector<Vec3f> read(grid.size());
    for (std::size_t i = 0; i < grid.size(); ++i) {
        read[i] = work.points.position(grid[i]);
    }
    std::vector<Vec3f> write = read;
    auto at = [columns](int r, int c) { return static_cast<std::size_t>(r * columns + c); };
    for (int iteration = 0; iteration < gridIterations; ++iteration) {
        for (int r = 1; r + 1 < rows; ++r) {
            for (int c = 1; c + 1 < columns; ++c) {
                const Vec3f along = read[at(r, c + 1)] * rowWeight + read[at(r, c - 1)] * (1.0f - rowWeight);
                const Vec3f across =
                    read[at(r + 1, c)] * columnWeight + read[at(r - 1, c)] * (1.0f - columnWeight);
                write[at(r, c)] = ((along + across) * 0.5f) * columnWeight + write[at(r, c)] * 0.5f;
            }
        }
        std::swap(read, write);
    }
    const int base = work.points.base();
    for (std::size_t i = 0; i < grid.size(); ++i) {
        if (grid[i] >= base) {
            work.points.move(grid[i], read[i]);
        }
    }
}

} // namespace areko::remesher::extract
