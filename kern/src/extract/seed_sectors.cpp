#include "seed_sectors.h"

namespace areko::remesher::extract {

using mesh::isWholeNumber;

namespace {

// Corner of `point` in triangle t, else the corner at `point` in a triangle
// across an edge of t (edges of t in corner order); -1 when there is none.
int startCorner(const TraceFrame& frame, int t, int point)
{
    if (const int c = frame.cornerAt(t, point); c >= 0) {
        return c;
    }
    const mesh::Adjacency& mesh = frame.mesh();
    for (int c = 3 * t; c < 3 * t + 3; ++c) {
        const int e = mesh.cornerEdge(c);
        const int other = e >= 0 ? mesh.otherCorner(e, c) : -1;
        if (other >= 0 && mesh.vertex(mesh.prev(other)) == point) {
            return mesh.prev(other);
        }
    }
    return -1;
}

} // namespace

bool SeedSectors::choose(const TraceFrame& frame, int point, int& triangle, int& direction, std::array<float, 2>& uv)
{
    const mesh::Adjacency& mesh = frame.mesh();
    const int first = triangle;
    int e = startCorner(frame, first, point);
    if (e < 0) {
        return false;
    }
    triangle = e / 3;
    const bool odd = (direction & 1) != 0;
    const int valence = static_cast<int>(mesh.pointCorners(point).size());
    int kept = -1;
    for (int steps = 0;; ++steps) {
        if (steps > valence + 2) {
            return false;
        }
        const std::array<float, 2>& p = frame.cornerUv(e);
        if (!isWholeNumber(p[0]) || !isWholeNumber(p[1])) {
            return false;
        }
        const std::array<float, 2>& n = frame.cornerUv(mesh.next(e));
        const std::array<float, 2>& r = frame.cornerUv(mesh.prev(e));
        const float tu = (p[0] - n[0]) / (r[0] - n[0]);
        const float tv = (p[1] - n[1]) / (r[1] - n[1]);
        const bool inU = tu >= 0.0f && tu <= 1.0f;
        const bool inV = tv >= 0.0f && tv <= 1.0f;
        const bool start = triangle == first;
        auto take = [&](int d) {
            direction = d;
            uv = p;
            return true;
        };
        if (!inU) {
            if (inV && (!start || odd)) {
                return take(0);
            }
        } else if (!start || !(tu == tv)) {
            if (!inV) {
                if (!start || !odd) {
                    return take(1);
                }
            } else if (tu <= tv) {
                if (!start) {
                    return take(0);
                }
                if (!odd) {
                    return take(1);
                }
                kept = 0;
            } else {
                if (!start) {
                    return take(1);
                }
                if (odd) {
                    return take(0);
                }
                kept = 1;
            }
        }
        // On to the triangle across the edge of corner e.
        const int edge = mesh.cornerEdge(e);
        const int other = edge >= 0 ? mesh.otherCorner(edge, e) : -1;
        if (other < 0) {
            triangle = -1;
            return true;
        }
        e = mesh.next(other);
        triangle = e / 3;
        if (triangle == first) {
            if (kept == -1) {
                direction = (direction + 1) % 4;
                return false;
            }
            direction = kept;
            return true;
        }
    }
}

} // namespace areko::remesher::extract
