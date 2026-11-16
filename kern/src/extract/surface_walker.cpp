#include "surface_walker.h"

#include "trace_node.h"

#include "field/axis_rotation.h"
#include "mesh/measures.h"

#include <cmath>

namespace areko::remesher::extract {

using mesh::Vec3d;
using mesh::Vec3f;

namespace {

// The turn of a direction from one triangle plane into the next: built and
// applied in double, from an axis that is not normalized first.
using Turn = field::AxisRotation<field::PlaneTurn>;

constexpr float startTriangle = -1e-6f;
constexpr float startFan = -100.0f;
constexpr float shortestStep = 1e-5f;
constexpr double smallestTurn = 0.001;
constexpr int maxSteps = 300;

Vec3f narrow(const Vec3d& v)
{
    return {static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z)};
}

// v without its part along n, in double.
void dropNormal(Vec3d& v, const Vec3d& n)
{
    const double t = mesh::dot(v, n);
    v = {v.x - n.x * t, v.y - n.y * t, v.z - n.z * t};
}

// Float vector against a double direction: products and sums in double, then float.
float along(Vec3f u, const Vec3d& d)
{
    return static_cast<float>((static_cast<double>(u.x) * d.x + static_cast<double>(u.y) * d.y) +
                              static_cast<double>(u.z) * d.z);
}

// a + (b - a) s in float.
Vec3f between(Vec3f a, Vec3f b, float s) { return a + (b - a) * s; }

} // namespace

SurfaceWalker::SurfaceWalker(const mesh::Mesh& triangles, const mesh::BitArray& fixedEdges)
    : adjacency_(triangles.adjacency()), points_(triangles.points()), normals_(triangles.triangleGeometry().normals()),
      fixed_(fixedEdges), useFixed_(fixedEdges.size() == static_cast<std::size_t>(adjacency_.edgeCount()))
{
}

Vec3f SurfaceWalker::normalOver(std::span<const int> corners) const
{
    Vec3f sum{};
    for (const int c : corners) {
        sum = sum + normals_[static_cast<std::size_t>(c / 3)];
    }
    mesh::normalizeInPlace(sum);
    return sum;
}

std::optional<SurfaceWalker::Crossing> SurfaceWalker::ray(Vec3f a, Vec3f b, const Vec3d& o, const Vec3d& d)
{
    const Vec3d e = mesh::wide(b - a);
    const double alpha = mesh::dot(e, d);
    const double beta = mesh::dot(e, e);
    const double den = beta - alpha * alpha;
    if (1e-50 > den) {
        return std::nullopt;
    }
    const Vec3d w = mesh::wide(a) - o;
    const double f6 = mesh::dot(d, w);
    const double f3 = mesh::dot(e, w);
    return Crossing{(beta * f6 - alpha * f3) / den, (f6 * alpha - f3) / den};
}

Vec3f SurfaceWalker::turned(Vec3f n1, Vec3f n2, Vec3f d)
{
    const Vec3f axis = mesh::cross(n1, n2);
    const float angle = mesh::vectorAngle(n1, n2);
    Turn turn;
    if (1e-7 > mesh::wide(mesh::dot(axis, axis))) {
        if (static_cast<double>(angle) > mesh::halfPiD) {
            turn = Turn(mesh::cross(mesh::smallestAxis(n1), n1), mesh::piF);
        }
    } else {
        turn = Turn(axis, angle);
    }
    // Rows summed in double, + 0, then rounded. The projective divisor
    // ((0 x + 0 y) + 0 z) + 1 is 1 for every finite d; a NaN in d makes every
    // row NaN either way, so the division is left out.
    return turn.apply(d);
}

void SurfaceWalker::walk(PointSource& source, Vec3f& position, Vec3f shift)
{
    if (shift.x == 0.0f && shift.y == 0.0f && shift.z == 0.0f) {
        return;
    }
    Vec3f normal;
    if (source.type == node_kind::point) {
        normal = normalOver(adjacency_.pointCorners(source.id));
    } else if (source.type == node_kind::edge) {
        normal = normalOver(adjacency_.edgeCorners(source.id));
    } else if (source.type == 4) {
        normal = normals_[static_cast<std::size_t>(source.id)];
    } else {
        return;
    }
    Vec3d d = mesh::wide(shift);
    dropNormal(d, mesh::wide(normal));
    double rest = mesh::normalizeInPlace(d);
    if (rest == 0.0) {
        return;
    }
    Vec3d o = mesh::wide(position);
    int entry = -1;
    for (int step = 1; step <= maxSteps; ++step) {
        Vec3f exit{};
        int exitCorner = -1;
        // Where the walk goes on after a stretch along a fixed edge or to a fixed edge's end.
        PointSource next{};
        int nextEntry = -1;
        bool fan = false;
        if (source.type == 4) {
            const int t = source.id;
            float best = startTriangle;
            for (int k = 0; k < 3; ++k) {
                const int c = 3 * t + k;
                if (adjacency_.cornerEdge(c) == entry) {
                    continue;
                }
                const Vec3f a = points_[static_cast<std::size_t>(adjacency_.vertex(c))];
                const Vec3f b = points_[static_cast<std::size_t>(adjacency_.vertex(3 * t + (k + 1) % 3))];
                const std::optional<Crossing> x = ray(a, b, o, d);
                if (x && 0.0 <= x->s && x->s <= 1.0 && x->t >= static_cast<double>(best)) {
                    const float s = static_cast<float>(x->s);
                    best = static_cast<float>(x->t);
                    exit = between(a, b, s);
                    exitCorner = c;
                }
            }
            if (exitCorner < 0) {
                position = narrow(o);
                return;
            }
            if (!(0.0f < best)) {
                // Already on the border of the triangle: its corner or side becomes the source.
                bool corner = false;
                for (int k = 0; k < 3 && !corner; ++k) {
                    const int v = adjacency_.vertex(3 * t + k);
                    if (exit == points_[static_cast<std::size_t>(v)]) {
                        source = {node_kind::point, v};
                        corner = true;
                    }
                }
                if (!corner) {
                    source = {node_kind::edge, adjacency_.cornerEdge(exitCorner)};
                }
                continue;
            }
            const double dist = static_cast<double>(best);
            if (rest <= dist) {
                position = {static_cast<float>(d.x * rest + o.x), static_cast<float>(d.y * rest + o.y),
                            static_cast<float>(d.z * rest + o.z)};
                return;
            }
            rest -= best < shortestStep ? static_cast<double>(shortestStep) : dist;
            o = mesh::wide(exit);
            const int e = adjacency_.cornerEdge(exitCorner);
            entry = e;
            const int across = e >= 0 ? adjacency_.otherCorner(e, exitCorner) : -1;
            if (across < 0) {
                position = exit;
                return;
            }
            const int beyond = across / 3;
            source.id = beyond;
            Vec3d laid = d;
            dropNormal(laid, mesh::wide(normals_[static_cast<std::size_t>(beyond)]));
            if (smallestTurn > mesh::normalizeInPlace(laid)) {
                laid = mesh::wide(turned(normals_[static_cast<std::size_t>(exitCorner / 3)],
                                         normals_[static_cast<std::size_t>(beyond)], narrow(d)));
            }
            if (rest <= 0.0) {
                break;
            }
            d = laid;
            continue;
        }
        if (source.type == node_kind::point) {
            const int v = source.id;
            bool alongFixed = false;
            if (useFixed_) {
                int count = 0;
                for (const int e : adjacency_.pointEdges(v)) {
                    count += fixed(e) ? 1 : 0;
                }
                if (count != 0 && count != 2) {
                    position = narrow(o);
                    return;
                }
                if (count == 2) {
                    // Along the fixed edge that points most along the walk, to its other end.
                    float bestAlong = -1.0f;
                    int chosen = -1;
                    for (const int e : adjacency_.pointEdges(v)) {
                        if (!fixed(e)) {
                            continue;
                        }
                        const Vec3f end = points_[static_cast<std::size_t>(adjacency_.otherEnd(e, v))];
                        Vec3f r = end - narrow(o);
                        mesh::normalizeInPlace(r);
                        const float a = along(r, d);
                        if (bestAlong < a) {
                            bestAlong = a;
                            exit = end;
                            chosen = e;
                        }
                    }
                    if (chosen < 0) {
                        position = narrow(o);
                        return;
                    }
                    next = {node_kind::edge, chosen};
                    alongFixed = true;
                }
            }
            if (!alongFixed) {
                // Out through the side opposite v of one triangle of the fan.
                const auto corners = adjacency_.pointCorners(v);
                if (corners.empty()) {
                    break;
                }
                float best = startFan;
                for (const int c : corners) {
                    const int cn = adjacency_.next(c);
                    const Vec3f a = points_[static_cast<std::size_t>(adjacency_.vertex(cn))];
                    const Vec3f b = points_[static_cast<std::size_t>(adjacency_.vertex(adjacency_.prev(c)))];
                    const std::optional<Crossing> x = ray(a, b, o, d);
                    if (x && 0.0 <= x->s && x->s <= 1.0 && x->t > 0.0 && x->t > static_cast<double>(best)) {
                        best = static_cast<float>(x->t);
                        exit = between(a, b, static_cast<float>(x->s));
                        exitCorner = cn;
                    }
                }
                if (exitCorner < 0) {
                    break;
                }
                fan = true;
            }
        } else if (source.type == node_kind::edge) {
            const int e = source.id;
            if (fixed(e)) {
                // Only along the edge, to the end it points to; R shrinks with the share along it.
                const auto ends = adjacency_.edgeVertices(e);
                Vec3f u = points_[static_cast<std::size_t>(ends[1])] - points_[static_cast<std::size_t>(ends[0])];
                mesh::normalizeInPlace(u);
                const float a = along(u, d);
                rest = mesh::wide(std::fabs(a)) * rest;
                const int w = a < 0.0f ? ends[0] : ends[1];
                if (w < 0) {
                    break;
                }
                exit = points_[static_cast<std::size_t>(w)];
                next = {node_kind::point, w};
                nextEntry = e;
            } else {
                // Out through one of the two other sides of a triangle at the edge.
                const auto corners = adjacency_.edgeCorners(e);
                if (corners.empty()) {
                    break;
                }
                float best = startFan;
                auto test = [&](int from, int to, int corner) {
                    const Vec3f a = points_[static_cast<std::size_t>(adjacency_.vertex(from))];
                    const Vec3f b = points_[static_cast<std::size_t>(adjacency_.vertex(to))];
                    const std::optional<Crossing> x = ray(a, b, o, d);
                    if (x && 0.0 <= x->s && x->s <= 1.0 && x->t > 0.0 && x->t > static_cast<double>(best)) {
                        best = static_cast<float>(x->t);
                        exit = between(a, b, static_cast<float>(x->s));
                        exitCorner = corner;
                    }
                };
                for (const int c : corners) {
                    const int cn = adjacency_.next(c);
                    const int cp = adjacency_.prev(c);
                    test(cn, cp, cn);
                    test(cp, c, cp);
                }
                if (exitCorner < 0) {
                    break;
                }
                fan = true;
            }
        } else {
            break;
        }

        // The stretch to `exit`.
        const Vec3d gap = o - mesh::wide(exit);
        const float dist = static_cast<float>(mesh::length(gap));
        if (rest <= static_cast<double>(dist)) {
            const Vec3f from = narrow(o);
            if (dist == 0.0f) {
                position = from;
                return;
            }
            const float share = static_cast<float>(rest / static_cast<double>(dist));
            position = (exit - from) * share + from;
            return;
        }
        rest -= dist < shortestStep ? static_cast<double>(shortestStep) : static_cast<double>(dist);
        o = mesh::wide(exit);
        if (fan) {
            const int e = adjacency_.cornerEdge(exitCorner);
            entry = e;
            const int across = e >= 0 ? adjacency_.otherCorner(e, exitCorner) : -1;
            if (across < 0) {
                position = exit;
                return;
            }
            source = {4, across / 3};
            dropNormal(d, mesh::wide(normals_[static_cast<std::size_t>(across / 3)]));
            (void)mesh::normalizeInPlace(d);
        } else {
            source = next;
            entry = nextEntry;
        }
        if (rest <= 0.0) {
            break;
        }
    }
    position = narrow(o);
}

} // namespace areko::remesher::extract
