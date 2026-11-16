#include "decimator.h"

#include "mesh/sharp_edges.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr float maxCost = std::numeric_limits<float>::max();
constexpr float lockedCost = std::bit_cast<float>(0x5368d4a5u); // 999999995904
constexpr double lockedCostD = 999999995904.0;
constexpr float flipDotEps = std::bit_cast<float>(0x26901d7du);   // 1e-15
constexpr float flipCos2Limit = std::bit_cast<float>(0x3bf8e8e8u); // cos^2(85 degrees)
constexpr float detEps = std::bit_cast<float>(0x3727c5acu);       // 1e-5
constexpr double lineAngle = 0.3490658503988659;                   // 20 degrees
constexpr int minVertexCount = 50;
constexpr int progressInterval = 5001;
// Blend factors of the fallback positions between the two vertices.
constexpr std::array<float, 9> fallbackSteps = {0.0f, 1.0f, 0.5f, 0.125f, 0.25f, 0.375f, 0.625f, 0.75f, 0.875f};

constexpr std::uint8_t activeBit = 1;
constexpr std::uint8_t borderBit = 2;
constexpr std::uint8_t lockedBit = 4;
constexpr std::uint8_t lineBits = 0xc;
constexpr std::uint8_t triangleActive = 8;

std::uint8_t channelByte(float v)
{
    return static_cast<std::uint8_t>(truncToInt(v));
}

} // namespace

Decimator::Decimator(float lengthScale, int targetVertexCount, bool useSharpEdges, run::Progress& progress)
    : target_(targetVertexCount), useSharpEdges_(useSharpEdges), progress_(progress)
{
    const float l = lengthScale / 6.0f;
    const float l2 = l * l;
    penaltyScale_ = detEps;
    lengthLimit_ = static_cast<float>(static_cast<double>(l2) * 1e-8);
    shortEdge2_ = static_cast<float>(static_cast<double>(l2) * 1e-6);
    penaltyOffset_ = shortEdge2_ / 3.0f;
}

Mesh Decimator::run(Mesh&& source)
{
    load(source);
    source = Mesh();
    computeFlags();
    startCosts();
    collapseLoop();
    return output();
}

void Decimator::load(Mesh& source)
{
    tri_ = source.topology().cornerVertex();
    pos_ = source.points();
    const std::size_t vertices = pos_.size();
    const std::size_t triangles = tri_.size() / 3;
    AttributeMaps& maps = source.maps();
    if (const auto* c4 = maps.findAs<TypedMap<Color4ub>>(map_names::colors, ElementKind::Point)) {
        colors_ = c4->values();
    } else if (const auto* c3 = maps.findAs<TypedMap<Color3ub>>(map_names::colors, ElementKind::Point)) {
        for (const Color3ub& c : c3->values()) {
            colors_.push_back({c.r, c.g, c.b, 255});
        }
    }
    if (colors_.size() != vertices) {
        colors_.clear();
    }
    mixColors_ = !colors_.empty();
    if (const auto* d = maps.findAs<TypedMap<float>>(map_names::density, ElementKind::Point);
        d != nullptr && d->size() == vertices) {
        density_ = d->values();
    }
    if (const auto* g = maps.findAs<TypedMap<int>>(map_names::faceGroups, ElementKind::Face); g != nullptr && g->size() == triangles) {
        groups_ = g->values();
        hasGroups_ = true;
    }
    triFlags_.assign(triangles, triangleActive);
    if (useSharpEdges_) {
        SharpEdges::toHalfEdgeForm(source);
        if (const auto* sharp = maps.findAs<BoolMap>(map_names::sharpEdges, ElementKind::HalfEdge);
            sharp != nullptr && sharp->bits().size() == 3 * triangles) {
            for (std::size_t t = 0; t < triangles; ++t) {
                const BitArray& b = sharp->bits();
                triFlags_[t] = static_cast<std::uint8_t>(triangleActive | (b.test(3 * t) ? 1 : 0) | (b.test(3 * t + 1) ? 2 : 0) |
                                                         (b.test(3 * t + 2) ? 4 : 0));
            }
            hasSharp_ = true;
        }
    }
    Vec3f lo{};
    Vec3f hi{};
    for (std::size_t v = 0; v < vertices; ++v) {
        const Vec3f p = pos_[v];
        if (v == 0) {
            lo = p;
            hi = p;
            continue;
        }
        lo = {p.x < lo.x ? p.x : lo.x, p.y < lo.y ? p.y : lo.y, p.z < lo.z ? p.z : lo.z};
        hi = {p.x > hi.x ? p.x : hi.x, p.y > hi.y ? p.y : hi.y, p.z > hi.z ? p.z : hi.z};
    }
    const Vec3f d = hi - lo;
    diagonal_ = std::sqrt((d.x * d.x + d.y * d.y) + d.z * d.z);
}

std::vector<int> Decimator::cornersOf(int v) const
{
    std::vector<int> out;
    for (int c = head_[static_cast<std::size_t>(v)]; c >= 0; c = next_[static_cast<std::size_t>(c)]) {
        out.push_back(c);
    }
    return out;
}

std::vector<int> Decimator::neighbours(int v) const
{
    std::vector<int> out;
    auto add = [&out](int n) {
        if (std::find(out.begin(), out.end(), n) == out.end()) {
            out.push_back(n);
        }
    };
    for (int c = head_[static_cast<std::size_t>(v)]; c >= 0; c = next_[static_cast<std::size_t>(c)]) {
        add(tri_[static_cast<std::size_t>(nextCorner(c))]);
        add(tri_[static_cast<std::size_t>(prevCorner(c))]);
    }
    return out;
}

bool Decimator::sharpBit(int halfEdge) const noexcept
{
    return ((triFlags_[static_cast<std::size_t>(halfEdge / 3)] >> (halfEdge % 3)) & 1) != 0;
}

void Decimator::unlink(int vertex, int corner)
{
    int* link = &head_[static_cast<std::size_t>(vertex)];
    while (*link >= 0) {
        if (*link == corner) {
            *link = next_[static_cast<std::size_t>(corner)];
            return;
        }
        link = &next_[static_cast<std::size_t>(*link)];
    }
}

void Decimator::pushFront(int vertex, int corner)
{
    next_[static_cast<std::size_t>(corner)] = head_[static_cast<std::size_t>(vertex)];
    head_[static_cast<std::size_t>(vertex)] = corner;
}

void Decimator::computeFlags()
{
    const std::size_t vertices = pos_.size();
    head_.assign(vertices, -1);
    next_.assign(tri_.size(), -1);
    for (std::size_t c = 0; c < tri_.size(); ++c) {
        pushFront(tri_[c], static_cast<int>(c));
    }
    flags_.assign(vertices, 0);
    int excluded = 0;
    struct Pair {
        int vertex;
        int first;
        int second;
    };
    std::vector<Pair> pairs;
    std::vector<Vec3f> directions;
    for (std::size_t i = 0; i < vertices; ++i) {
        pairs.clear();
        int corners = 0;
        auto add = [&pairs](int n, int id) {
            for (Pair& p : pairs) {
                if (p.vertex == n) {
                    p.second = id;
                    return;
                }
            }
            pairs.push_back({n, id, -1});
        };
        for (int c = head_[i]; c >= 0; c = next_[static_cast<std::size_t>(c)]) {
            ++corners;
            add(tri_[static_cast<std::size_t>(nextCorner(c))], c);
            add(tri_[static_cast<std::size_t>(prevCorner(c))], prevCorner(c));
        }
        if (corners == 0) {
            ++excluded;
            continue;
        }
        if (corners != static_cast<int>(pairs.size())) {
            flags_[i] = 0xf;
            continue;
        }
        std::uint8_t flag = activeBit;
        if (hasGroups_ || hasSharp_) {
            directions.clear();
            for (const Pair& p : pairs) {
                if (p.second < 0) {
                    continue;
                }
                const bool marked =
                    (hasGroups_ && groups_[static_cast<std::size_t>(p.first / 3)] != groups_[static_cast<std::size_t>(p.second / 3)]) ||
                    (hasSharp_ && (sharpBit(p.first) || sharpBit(p.second)));
                if (marked) {
                    directions.push_back(pos_[static_cast<std::size_t>(p.vertex)] - pos_[i]);
                }
            }
            if (directions.size() >= 3) {
                flag |= lineBits;
            } else if (directions.size() == 2) {
                const Vec3f d1 = directions[0];
                const Vec3f d2 = directions[1];
                const double l1 = lengthAsDouble(d1);
                const double l2 = lengthAsDouble(d2);
                const double lengths = l2 * l1;
                bool corner = false;
                if (lengths != 0.0) {
                    const double cosine = static_cast<double>(dot(d1, d2)) / lengths;
                    // A straight line runs through the point: the two directions point away from each other.
                    const double c = -cosine;
                    if (!(c >= 1.0)) {
                        const double a = static_cast<double>(static_cast<float>(std::acos(c)));
                        corner = !(a <= lineAngle);
                    }
                }
                flag |= corner ? lineBits : std::uint8_t{8};
            }
        }
        flags_[i] = flag;
    }
    active_ = static_cast<int>(vertices) - excluded;
}

Quadric Decimator::triangleQuadric(int triangle) const
{
    const auto t = static_cast<std::size_t>(triangle) * 3;
    return Quadric::ofTriangle(pos_[static_cast<std::size_t>(tri_[t])], pos_[static_cast<std::size_t>(tri_[t + 1])],
                               pos_[static_cast<std::size_t>(tri_[t + 2])], weightScale_);
}

Quadric Decimator::vertexQuadric(const std::vector<int>& corners) const
{
    Quadric q;
    for (const int c : corners) {
        q += triangleQuadric(c / 3);
    }
    return q;
}

void Decimator::fallback(int a, int b, const Quadric& q, bool lineCase, bool chooseB, float detAbs, Evaluation& e) const
{
    const Vec3f pa = pos_[static_cast<std::size_t>(a)];
    const Vec3f pb = pos_[static_cast<std::size_t>(b)];
    if (lineCase) {
        e.position = chooseB ? pb : pa;
        e.cost = static_cast<float>(q.evaluate(e.position));
        return;
    }
    float best = e.cost;
    for (std::size_t i = 0; i < fallbackSteps.size(); ++i) {
        const float t = fallbackSteps[i];
        const float s = 1.0f - t;
        const Vec3f p{pb.x * t + pa.x * s, pb.y * t + pa.y * s, pb.z * t + pa.z * s};
        const float c = static_cast<float>(q.evaluate(p));
        if (i == 2) {
            if (static_cast<double>(detAbs) < 1e-10 && static_cast<double>(best) >= static_cast<double>(c) * 0.999) {
                e.cost = c;
                e.position = p;
                best = c;
                if (static_cast<double>(best) < 1e-15) {
                    break;
                }
                continue;
            }
            if (static_cast<double>(best) < 1e-15) {
                break;
            }
        }
        if (best > c) {
            e.cost = c;
            e.position = p;
            best = c;
        }
    }
    if (0.0f > best) {
        e.cost = 0.0f;
    }
}

Decimator::Evaluation Decimator::evaluate(int a, int b, const Quadric& quadric) const
{
    Evaluation e;
    e.position = pos_[static_cast<std::size_t>(a)];
    const std::uint8_t fa = flags_[static_cast<std::size_t>(a)];
    const std::uint8_t fb = flags_[static_cast<std::size_t>(b)];
    const std::uint8_t f = fa | fb;
    if (f & lockedBit) {
        e.cost = maxCost;
        return e;
    }
    if (f & lineBits) {
        e.cost = maxCost;
        fallback(a, b, quadric, true, (fa & lineBits) == 0, 0.0f, e);
    } else {
        const auto& q = quadric.q;
        const double c00 = q[3] * q[5] - q[4] * q[4];
        const double m01 = q[1] * q[5] - q[4] * q[2];
        const double m02 = q[4] * q[1] - q[3] * q[2];
        const double det = (q[0] * c00 - q[1] * m01) + q[2] * m02;
        if (det == 0.0) {
            e.cost = maxCost;
            fallback(a, b, quadric, false, false, std::fabs(static_cast<float>(det)), e);
        } else {
            const double i00 = c00 / det;
            const double i01 = (-m01) / det;
            const double i02 = m02 / det;
            const double i11 = (q[0] * q[5] - q[2] * q[2]) / det;
            const double i12 = (-(q[0] * q[4] - q[2] * q[1])) / det;
            const double i22 = (q[0] * q[3] - q[1] * q[1]) / det;
            e.position = {static_cast<float>(-((i00 * q[6] + i01 * q[7]) + i02 * q[8])),
                          static_cast<float>(-((i01 * q[6] + i11 * q[7]) + i12 * q[8])),
                          static_cast<float>(-((i02 * q[6] + i12 * q[7]) + i22 * q[8]))};
            const float c = static_cast<float>(quadric.evaluate(e.position));
            e.cost = 0.0f > c ? 0.0f : c;
            const float detAbs = std::fabs(static_cast<float>(det));
            if (detEps > detAbs) {
                fallback(a, b, quadric, false, false, detAbs, e);
            }
        }
    }
    const float c = e.cost;
    if (lengthLimit_ > c) {
        const Vec3f d = pos_[static_cast<std::size_t>(a)] - pos_[static_cast<std::size_t>(b)];
        const float l2 = (d.x * d.x + d.y * d.y) + d.z * d.z;
        if (l2 > shortEdge2_) {
            e.cost = (l2 - penaltyOffset_) * penaltyScale_ + c;
        }
    }
    return e;
}

void Decimator::startCosts()
{
    const std::size_t vertices = pos_.size();
    const std::size_t triangles = tri_.size() / 3;
    weightScale_ = static_cast<float>(
        1.0 / ((static_cast<double>(diagonal_) * static_cast<double>(diagonal_)) / static_cast<double>(triangles)));
    std::vector<Quadric> quadrics(vertices);
    for (std::size_t v = 0; v < vertices; ++v) {
        quadrics[v] = vertexQuadric(cornersOf(static_cast<int>(v)));
    }
    best_.assign(vertices, Collapse{-1, maxCost});
    for (std::size_t v = 0; v < vertices; ++v) {
        Collapse& best = best_[v];
        for (const int n : neighbours(static_cast<int>(v))) {
            const Evaluation e = evaluate(static_cast<int>(v), n, quadrics[v] + quadrics[static_cast<std::size_t>(n)]);
            if (e.cost < best.cost) {
                best = {n, e.cost};
            }
        }
    }
    queue_.emplace(best_, diagonal_);
    localOf_.assign(vertices, -1);
    triangleMark_.assign(triangles, 0);
}

int Decimator::check(int v, int t, Vec3f position, const std::vector<int>& cornersV, const std::vector<int>& cornersT) const
{
    std::vector<int> moved;
    for (const int c : cornersV) {
        moved.push_back(c / 3);
    }
    for (const int c : cornersT) {
        const int s = c / 3;
        const auto it = std::find(moved.begin(), moved.end(), s);
        if (it != moved.end()) {
            moved.erase(it);
        } else {
            moved.push_back(s);
        }
    }
    for (const int s : moved) {
        const auto b = static_cast<std::size_t>(s) * 3;
        std::array<Vec3f, 3> p;
        for (std::size_t j = 0; j < 3; ++j) {
            p[j] = pos_[static_cast<std::size_t>(tri_[b + j])];
        }
        const Vec3f e1 = p[1] - p[0];
        const Vec3f e2 = p[2] - p[0];
        const Vec3f o = cross(e1, e2);
        for (std::size_t j = 0; j < 3; ++j) {
            const int x = tri_[b + j];
            if (x == v || x == t) {
                p[j] = position;
            }
        }
        const Vec3f f1 = p[1] - p[0];
        const Vec3f f2 = p[2] - p[0];
        const Vec3f n = cross(f1, f2);
        const float d = (o.x * n.x + o.y * n.y) + o.z * n.z;
        if (d > flipDotEps) {
            continue;
        }
        // Also an undefined (NaN) orientation counts as folded over.
        if (!(d > 0.0f)) {
            return 3;
        }
        const float q = (d * d) / (((o.x * o.x + o.y * o.y) + o.z * o.z) * ((n.x * n.x + n.y * n.y) + n.z * n.z));
        if (flipCos2Limit > q) {
            return 3;
        }
    }
    return 0;
}

int Decimator::collapse(int v, int t, Vec3f position, std::vector<int> cornersV, Record& record,
                        std::vector<int>& ringCorners, std::vector<int>& ring)
{
    const std::uint8_t ft = flags_[static_cast<std::size_t>(t)];
    flags_[static_cast<std::size_t>(v)] &= static_cast<std::uint8_t>(~activeBit);
    --active_;
    const std::uint8_t fv = flags_[static_cast<std::size_t>(v)];
    record.removedTriangles.clear();
    bool nonManifold = false;
    if (!cornersV.empty()) {
        next_[static_cast<std::size_t>(cornersV.back())] = head_[static_cast<std::size_t>(t)];
        head_[static_cast<std::size_t>(t)] = head_[static_cast<std::size_t>(v)];
        head_[static_cast<std::size_t>(v)] = -1;
        for (const int c : cornersV) {
            tri_[static_cast<std::size_t>(c)] = t;
            const int a = tri_[static_cast<std::size_t>(nextCorner(c))];
            const int b = tri_[static_cast<std::size_t>(prevCorner(c))];
            if (a == t || b == t) {
                const int s = c / 3;
                triFlags_[static_cast<std::size_t>(s)] &= static_cast<std::uint8_t>(~triangleActive);
                record.removedTriangles.push_back(s);
                unlink(t, c);
                unlink(a, nextCorner(c));
                unlink(b, prevCorner(c));
                if ((head_[static_cast<std::size_t>(a)] < 0 && a != t) || (head_[static_cast<std::size_t>(b)] < 0 && b != t)) {
                    nonManifold = true;
                }
            }
        }
    }
    record.vertex = v;
    record.oldPosition = pos_[static_cast<std::size_t>(t)];
    record.corners = std::move(cornersV);

    int result = 0;
    if (nonManifold) {
        result = 3;
    } else {
        ringCorners = cornersOf(t);
        ring.clear();
        auto add = [&ring, t](int n) {
            if (n != t && std::find(ring.begin(), ring.end(), n) == ring.end()) {
                ring.push_back(n);
            }
        };
        for (const int c : ringCorners) {
            add(tri_[static_cast<std::size_t>(nextCorner(c))]);
            add(tri_[static_cast<std::size_t>(prevCorner(c))]);
        }
        const auto& removed = record.removedTriangles;
        if ((ft & 8) && (fv & 8) && hasGroups_) {
            if (removed.size() <= 1) {
                result = 3;
            } else {
                const int g = groups_[static_cast<std::size_t>(removed[0])];
                const bool same = std::all_of(removed.begin(), removed.end(),
                                              [&](int s) { return groups_[static_cast<std::size_t>(s)] == g; });
                result = same ? 3 : 0;
            }
        }
        if (result == 0) {
            const bool tBorder = (ft & borderBit) != 0;
            const bool vBorder = (fv & borderBit) != 0;
            const std::size_t a = ringCorners.size();
            const std::size_t b = ring.size();
            if (!tBorder && !vBorder) {
                result = a > b ? 3 : 0;
            } else if (!tBorder) {
                result = a >= b ? 3 : 0;
                flags_[static_cast<std::size_t>(t)] |= borderBit;
            } else if (!vBorder) {
                result = a >= b ? 3 : 0;
            } else {
                result = a + 1 != b ? 3 : 0;
                flags_[static_cast<std::size_t>(t)] |= borderBit;
            }
        }
    }
    flags_[static_cast<std::size_t>(t)] = static_cast<std::uint8_t>((flags_[static_cast<std::size_t>(v)] & lineBits) |
                                                                    flags_[static_cast<std::size_t>(t)]);
    pos_[static_cast<std::size_t>(t)] = position;
    if (result != 0) {
        undo(t, record, ft);
        return result;
    }
    if (hasSharp_ && (flags_[static_cast<std::size_t>(t)] & lineBits)) {
        std::vector<int> partners;
        for (const int s : record.removedTriangles) {
            const auto b = static_cast<std::size_t>(s) * 3;
            for (std::size_t j = 0; j < 3; ++j) {
                const int p = tri_[b + j];
                const int q = tri_[b + (j + 1) % 3];
                if (((triFlags_[static_cast<std::size_t>(s)] >> j) & 1) != 0 && p != q) {
                    const int w = p == t ? q : p;
                    if (std::find(partners.begin(), partners.end(), w) == partners.end()) {
                        partners.push_back(w);
                    }
                    break;
                }
            }
        }
        if (!partners.empty()) {
            auto isPartner = [&](int x) { return std::find(partners.begin(), partners.end(), x) != partners.end(); };
            for (int c = head_[static_cast<std::size_t>(t)]; c >= 0; c = next_[static_cast<std::size_t>(c)]) {
                const auto s = static_cast<std::size_t>(c / 3);
                if (isPartner(tri_[static_cast<std::size_t>(nextCorner(c))])) {
                    triFlags_[s] |= static_cast<std::uint8_t>(1 << (c % 3));
                }
                if (isPartner(tri_[static_cast<std::size_t>(prevCorner(c))])) {
                    triFlags_[s] |= static_cast<std::uint8_t>(1 << ((c + 2) % 3));
                }
            }
        }
    }
    return 0;
}

void Decimator::undo(int t, Record& record, std::uint8_t oldTargetFlag)
{
    const int v = record.vertex;
    pos_[static_cast<std::size_t>(t)] = record.oldPosition;
    flags_[static_cast<std::size_t>(t)] = oldTargetFlag;
    for (const int c : record.corners) {
        const int s = c / 3;
        if ((triFlags_[static_cast<std::size_t>(s)] & triangleActive) == 0) {
            for (int j = 0; j < 3; ++j) {
                const int corner = 3 * s + j;
                const int x = tri_[static_cast<std::size_t>(corner)];
                flags_[static_cast<std::size_t>(x)] |= activeBit;
                pushFront(x, corner);
            }
            triFlags_[static_cast<std::size_t>(s)] |= triangleActive;
        }
        unlink(t, c);
        pushFront(v, c);
        tri_[static_cast<std::size_t>(c)] = v;
    }
    flags_[static_cast<std::size_t>(v)] |= activeBit;
    ++active_;
    record = Record();
}

void Decimator::refresh(int t, const std::vector<int>& ringCorners, const std::vector<int>& ring, std::vector<Update>& updates)
{
    std::vector<int> local{t};
    local.insert(local.end(), ring.begin(), ring.end());
    for (std::size_t k = 0; k < local.size(); ++k) {
        localOf_[static_cast<std::size_t>(local[k])] = static_cast<int>(k);
    }
    std::vector<int> triangles;
    auto take = [&](int c) {
        const auto s = static_cast<std::size_t>(c / 3);
        if (!triangleMark_[s]) {
            triangleMark_[s] = 1;
            triangles.push_back(c / 3);
        }
    };
    for (const int c : ringCorners) {
        take(c);
    }
    std::vector<std::vector<int>> ringNeighbours;
    for (const int r : ring) {
        ringNeighbours.push_back(neighbours(r));
        for (const int n : ringNeighbours.back()) {
            if (localOf_[static_cast<std::size_t>(n)] < 0) {
                localOf_[static_cast<std::size_t>(n)] = static_cast<int>(local.size());
                local.push_back(n);
                for (int c = head_[static_cast<std::size_t>(n)]; c >= 0; c = next_[static_cast<std::size_t>(c)]) {
                    take(c);
                }
            }
        }
    }
    // As in the procedure, the quadrics take only the triangles at the corners of the ring
    // and at the ring's outer neighbours; triangles spanned by ring vertices alone stay out
    // (taking them in made the decimated mesh differ and lost a singularity at one eye).
    std::vector<Quadric> quadrics(local.size());
    for (const int s : triangles) {
        const Quadric q = triangleQuadric(s);
        for (int j = 0; j < 3; ++j) {
            const int l = localOf_[static_cast<std::size_t>(tri_[static_cast<std::size_t>(3 * s + j)])];
            if (l >= 0) {
                quadrics[static_cast<std::size_t>(l)] += q;
            }
        }
    }
    const int ringEnd = static_cast<int>(ring.size()) + 1;
    std::vector<Update> records(static_cast<std::size_t>(ringEnd));
    for (int k = 0; k < ringEnd; ++k) {
        records[static_cast<std::size_t>(k)] = {local[static_cast<std::size_t>(k)], -1, maxCost, -1};
    }
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const int r = ring[i];
        const int u = static_cast<int>(i) + 1;
        for (const int n : ringNeighbours[i]) {
            const int w = localOf_[static_cast<std::size_t>(n)];
            if (!(w == 0 || w >= ringEnd || r <= n)) {
                continue;
            }
            const Evaluation e = evaluate(r, n, quadrics[static_cast<std::size_t>(u)] + quadrics[static_cast<std::size_t>(w)]);
            Update& ru = records[static_cast<std::size_t>(u)];
            if (e.cost < ru.cost) {
                ru = {r, n, e.cost, queue_->bin(e.cost)};
            }
            if (w < ringEnd) {
                Update& rw = records[static_cast<std::size_t>(w)];
                if (e.cost < rw.cost) {
                    rw = {n, r, e.cost, queue_->bin(e.cost)};
                }
            }
        }
    }
    for (const int x : local) {
        localOf_[static_cast<std::size_t>(x)] = -1;
    }
    for (const int s : triangles) {
        triangleMark_[static_cast<std::size_t>(s)] = 0;
    }
    updates.insert(updates.end(), records.begin(), records.end());
}

float Decimator::mixShare(int v, int t, Vec3f position) const
{
    const float dv = distance(pos_[static_cast<std::size_t>(v)], position);
    const float dt = distance(pos_[static_cast<std::size_t>(t)], position);
    return dv / (dt + dv);
}

float Decimator::mixDensity(int v, int t, float lambda) const
{
    const float mu = 1.0f - lambda;
    return density_[static_cast<std::size_t>(t)] * lambda + density_[static_cast<std::size_t>(v)] * mu;
}

Color4ub Decimator::mixColor(int v, int t, float lambda) const
{
    const Color4ub cv = colors_[static_cast<std::size_t>(v)];
    const Color4ub ct = colors_[static_cast<std::size_t>(t)];
    const float mu = 1.0f - lambda;
    auto rgb = [&](std::uint8_t a, std::uint8_t b) {
        return channelByte(((static_cast<float>(b) / 255.0f) * lambda + (static_cast<float>(a) / 255.0f) * mu) * 255.0f);
    };
    const std::uint8_t alpha =
        channelByte(((static_cast<float>(cv.a) * mu) / 255.0f + (static_cast<float>(ct.a) * lambda) / 255.0f) * 255.0f);
    return {rgb(cv.r, ct.r), rgb(cv.g, ct.g), rgb(cv.b, ct.b), alpha};
}

void Decimator::collapseLoop()
{
    const int startCount = active_;
    const float span = static_cast<float>(std::max(1, startCount - target_));
    const int lockedBin = queue_->binOfDouble(lockedCostD);
    std::vector<Update> updates;
    std::vector<int> ringCorners;
    std::vector<int> ring;
    Record record;
    int v = queue_->front();
    for (int round = 0;; ++round) {
        if (v < 0) {
            break;
        }
        const float k = best_[static_cast<std::size_t>(v)].cost;
        const int t = best_[static_cast<std::size_t>(v)].target;
        if (-1.0f < k || std::isnan(k)) {
            if (active_ <= target_) {
                break;
            }
        } else if (active_ <= minVertexCount) {
            break;
        }
        updates.clear();
        int removeOnly = -1;
        bool stop = false;
        auto lock = [&](int x) {
            if (queue_->queued(x)) {
                updates.push_back({x, best_[static_cast<std::size_t>(x)].target, lockedCost, lockedBin});
            } else {
                best_[static_cast<std::size_t>(x)].cost = lockedCost;
            }
        };
        if (t < 0 || (flags_[static_cast<std::size_t>(t)] & activeBit) == 0) {
            removeOnly = v;
        } else if (k >= lockedCost) {
            stop = true;
        } else {
            std::vector<int> cornersV = cornersOf(v);
            const std::vector<int> cornersT = cornersOf(t);
            const Vec3f position = evaluate(v, t, vertexQuadric(cornersT) + vertexQuadric(cornersV)).position;
            const bool withDensity = !density_.empty();
            const float lambda = mixColors_ || withDensity ? mixShare(v, t, position) : 0.0f;
            const Color4ub color = mixColors_ ? mixColor(v, t, lambda) : Color4ub{};
            const float density = withDensity ? mixDensity(v, t, lambda) : 0.0f;
            if (check(v, t, position, cornersV, cornersT) == 3) {
                lock(v);
            } else if (collapse(v, t, position, std::move(cornersV), record, ringCorners, ring) == 0) {
                removeOnly = v;
                if (mixColors_) {
                    colors_[static_cast<std::size_t>(t)] = color;
                }
                if (withDensity) {
                    density_[static_cast<std::size_t>(t)] = density;
                }
                refresh(t, ringCorners, ring, updates);
            } else {
                lock(v);
            }
        }
        for (const Update& u : updates) {
            queue_->remove(u.vertex);
            best_[static_cast<std::size_t>(u.vertex)] = {u.target, u.cost};
            if (u.cost != maxCost) {
                queue_->insert(u.vertex, u.bin);
            }
        }
        if (removeOnly >= 0) {
            queue_->remove(removeOnly);
        }
        v = queue_->front();
        if (round % progressInterval == 0) {
            progress_.check();
            progress_.report(static_cast<float>(startCount - active_) / span);
        }
        if (stop) {
            break;
        }
    }
}

Mesh Decimator::output() const
{
    const std::size_t triangles = tri_.size() / 3;
    std::vector<int> renumber(pos_.size(), -1);
    int count = 0;
    std::vector<int> corners;
    std::vector<int> groups;
    BitArray sharp;
    for (std::size_t s = 0; s < triangles; ++s) {
        if ((triFlags_[s] & triangleActive) == 0) {
            continue;
        }
        for (std::size_t j = 0; j < 3; ++j) {
            int& n = renumber[static_cast<std::size_t>(tri_[3 * s + j])];
            if (n < 0) {
                n = count++;
            }
            corners.push_back(n);
        }
        if (hasGroups_) {
            groups.push_back(groups_[s]);
        }
    }
    std::vector<Vec3f> points(static_cast<std::size_t>(count));
    std::vector<Color3ub> colors;
    if (!colors_.empty()) {
        colors.resize(static_cast<std::size_t>(count));
    }
    std::vector<float> density(density_.empty() ? std::size_t{0} : static_cast<std::size_t>(count));
    for (std::size_t v = 0; v < pos_.size(); ++v) {
        const int n = renumber[v];
        if (n < 0) {
            continue;
        }
        points[static_cast<std::size_t>(n)] = pos_[v];
        if (!colors.empty()) {
            colors[static_cast<std::size_t>(n)] = {colors_[v].r, colors_[v].g, colors_[v].b};
        }
        if (!density.empty()) {
            density[static_cast<std::size_t>(n)] = density_[v];
        }
    }
    if (hasSharp_ && triangles > 0) {
        sharp.resize(corners.size());
        std::size_t out = 0;
        for (std::size_t s = 0; s < triangles; ++s) {
            if ((triFlags_[s] & triangleActive) == 0) {
                continue;
            }
            for (std::size_t j = 0; j < 3; ++j, ++out) {
                if ((triFlags_[s] >> j) & 1) {
                    sharp.set(out);
                }
            }
        }
    }
    Mesh mesh;
    mesh.swapPoints(points);
    mesh.editTopology().swapCorners(corners, FaceMode::Triangles);
    mesh.topologyChanged();
    if (!colors.empty()) {
        mesh.maps().set(std::make_unique<TypedMap<Color3ub>>(std::string(map_names::colors), ElementKind::Point, std::move(colors)));
    }
    if (!groups.empty()) {
        mesh.maps().set(std::make_unique<TypedMap<int>>(std::string(map_names::faceGroups), ElementKind::Face, std::move(groups)));
    }
    if (hasSharp_ && triangles > 0) {
        mesh.maps().set(std::make_unique<BoolMap>(std::string(map_names::sharpEdges), ElementKind::HalfEdge, std::move(sharp)));
    }
    if (!density.empty()) {
        mesh.maps().set(std::make_unique<TypedMap<float>>(std::string(map_names::density), ElementKind::Point, std::move(density)));
    }
    return mesh;
}

} // namespace areko::remesher::prep
