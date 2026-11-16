#include "segment_cutter.h"

#include "curve_projector.h"
#include "segment_preparation.h"

#include "mesh/modification.h"
#include "prep/mesh_repair.h"
#include "prep/polygon_triangulator.h"
#include "util/introsort.h"
#include "util/rounding.h"

#include <CDT.h>

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <exception>
#include <map>
#include <tuple>
#include <utility>

namespace areko::remesher::curves {

using namespace areko::remesher::mesh;

namespace {

// One entry of a face: two positions on its border (corner k = k, side k at
// parameter t = k + t; -1 = inside the face) and their mesh points.
struct Entry {
    double tA = -1.0;
    double tB = -1.0;
    int vA = -1;
    int vB = -1;
};

using Polygon = std::vector<int>;

bool isWhole(double x) { return x == std::floor(x); }


// Points, sources and polygons of one cut.
class Cutter {
public:
    Cutter(Mesh& mesh, std::span<const CurveSegment> pieces, const SegmentCutter::Options& options)
        : mesh_(mesh), adjacency_(mesh.adjacency()), topology_(mesh.topology()), pieces_(pieces), options_(options),
          pointCount_(mesh.pointCount()), faceCount_(mesh.faceCount())
    {
    }

    SegmentCutter::Result run();

    [[nodiscard]] Vec3f position(int v) const
    {
        return v < pointCount_ ? mesh_.points()[static_cast<std::size_t>(v)]
                               : newPoints_[static_cast<std::size_t>(v - pointCount_)];
    }
    int addPoint(Vec3f p)
    {
        newPoints_.push_back(p);
        return pointCount_ + static_cast<int>(newPoints_.size()) - 1;
    }
    [[nodiscard]] int pointCount() const noexcept { return pointCount_; }

private:
    void toMeshPoints();
    [[nodiscard]] bool needsCut(const CurveSegment& piece, int i) const;
    [[nodiscard]] bool perimeter(float t, int e, int f, double& u) const;
    void countEntries();
    void writeEntries();
    void mergePoints();
    void rebuildFaces();
    [[nodiscard]] Mesh triangulatedFaces(std::vector<int>& origin) const;
    bool apply();
    void repairChains();
    void addSourceRow(std::initializer_list<std::pair<int, float>> row);
    void addFaceRow(int f);

    Mesh& mesh_;
    const Adjacency& adjacency_;
    const Topology& topology_;
    std::span<const CurveSegment> pieces_;
    SegmentCutter::Options options_;
    int pointCount_;
    int faceCount_;
    std::vector<Vec3f> newPoints_;
    Csr sources_;
    std::vector<float> weights_;
    std::vector<std::vector<int>> chains_;
    std::vector<std::vector<int>> numbers_;
    std::vector<std::vector<int>> meshPoints_; // per piece and entry
    std::vector<bool> closed_;
    BitArray inside_;
    std::vector<int> entryStarts_;
    std::vector<Entry> entries_;
    std::vector<int> moved_;
    std::vector<Vec3f> movedTo_;
    // New polygons, their origin and the faces they replace.
    std::vector<Polygon> polygons_;
    std::vector<int> origins_;
    BitArray removed_;
};

void Cutter::addSourceRow(std::initializer_list<std::pair<int, float>> row)
{
    for (const auto& [v, w] : row) {
        sources_.entries.push_back(v);
        weights_.push_back(w);
    }
    sources_.starts.push_back(static_cast<int>(sources_.entries.size()));
}

void Cutter::addFaceRow(int f)
{
    const int k = topology_.faceSize(f);
    const float w = topology_.mode() == FaceMode::Triangles ? 0.33333334f
                    : topology_.mode() == FaceMode::Quads ? 0.25f
                                                           : 1.0f / static_cast<float>(k);
    for (int c = topology_.faceBegin(f); c < topology_.faceEnd(f); ++c) {
        sources_.entries.push_back(topology_.vertex(c));
        weights_.push_back(w);
    }
    sources_.starts.push_back(static_cast<int>(sources_.entries.size()));
}

// Phase A.
void Cutter::toMeshPoints()
{
    const std::size_t count = pieces_.size();
    chains_.assign(count, {});
    numbers_.assign(count, {});
    meshPoints_.assign(count, {});
    closed_.assign(count, false);
    for (std::size_t s = 0; s < count; ++s) {
        const CurveSegment& piece = pieces_[s];
        const int n = piece.size();
        const bool closed = SegmentPreparation::isClosed(adjacency_, piece);
        closed_[s] = closed;
        std::vector<int>& kp = meshPoints_[s];
        if (n < 2 || (closed && piece.points.size() == 1)) {
            kp.assign(static_cast<std::size_t>(n), -1);
            continue;
        }
        kp.assign(static_cast<std::size_t>(n), 0);
        const auto faces = static_cast<int>(piece.faces.size());
        std::vector<int>& chain = chains_[s];
        for (int i = 0; i < n; ++i) {
            const Vec3f p = piece.points[static_cast<std::size_t>(i)];
            const int e = piece.edges[static_cast<std::size_t>(i)];
            int v = -1;
            if (e < 0) {
                if (i == n - 1 && closed) {
                    v = chain.front();
                } else {
                    const int fi = i == n - 1 ? faces - 1 : i;
                    const int f = fi >= 0 && fi < faces ? piece.faces[static_cast<std::size_t>(fi)] : -1;
                    if (f < 0) {
                        // No face: neither a point nor a source row.
                        kp[static_cast<std::size_t>(i)] = -1;
                        continue;
                    }
                    v = addPoint(p);
                    const int c0 = topology_.faceBegin(f);
                    if (topology_.faceSize(f) < 3) {
                        sources_.starts.push_back(static_cast<int>(sources_.entries.size()));
                    } else {
                        const int i0 = topology_.vertex(c0);
                        const int i1 = topology_.vertex(c0 + 1);
                        const int i2 = topology_.vertex(c0 + 2);
                        const Vec3f p0 = position(i0);
                        const Vec3f e1 = position(i1) - p0;
                        const Vec3f e2 = position(i2) - p0;
                        const Vec3f d = p - p0;
                        const float d11 = dot(e1, e1);
                        const float d22 = dot(e2, e2);
                        const float d12 = dot(e1, e2);
                        const float inv = 1.0f / (d11 * d22 - d12 * d12);
                        const float u = (d22 * dot(e1, d) - d12 * dot(e2, d)) * inv;
                        const float w = (d11 * dot(e2, d) - d12 * dot(e1, d)) * inv;
                        addSourceRow({{i0, (1.0f - u) - w}, {i1, u}, {i2, w}});
                    }
                }
            } else {
                const float t = piece.t[static_cast<std::size_t>(i)];
                const auto ends = adjacency_.edgeVertices(e);
                v = t == 0.0f ? ends[0] : (t == 1.0f ? ends[1] : -1);
                if (v >= 0) {
                    if (options_.moveCornerPoints) {
                        moved_.push_back(v);
                        movedTo_.push_back(p);
                    }
                } else if (i == n - 1 && closed) {
                    v = chain.front();
                } else {
                    v = addPoint(p);
                    addSourceRow({{ends[0], 1.0f - t}, {ends[1], t}});
                }
            }
            kp[static_cast<std::size_t>(i)] = v;
            if (chain.empty() || chain.back() != v) {
                chain.push_back(v);
                numbers_[s].push_back(i);
            }
        }
    }
    if (!moved_.empty()) {
        mesh_.movePoints(moved_, movedTo_);
    }
}

bool Cutter::needsCut(const CurveSegment& piece, int i) const
{
    if (i + 1 >= piece.size()) {
        return true;
    }
    const int e0 = piece.edges[static_cast<std::size_t>(i)];
    const int e1 = piece.edges[static_cast<std::size_t>(i) + 1];
    const int a = SegmentPreparation::cornerPoint(adjacency_, e0, piece.t[static_cast<std::size_t>(i)]);
    const int b = SegmentPreparation::cornerPoint(adjacency_, e1, piece.t[static_cast<std::size_t>(i) + 1]);
    if (topology_.mode() == FaceMode::Triangles) {
        // Nothing to cut when both ends are mesh points (they share an edge of the face).
        return !(a >= 0 && b >= 0);
    }
    if (a < 0 || e1 < 0) {
        return true;
    }
    if (a == b) {
        return false;
    }
    return adjacency_.edgeBetween(a, b) < 0;
}

bool Cutter::perimeter(float t, int e, int f, double& u) const
{
    const int c0 = topology_.faceBegin(f);
    const int k = topology_.faceSize(f);
    const int w = SegmentPreparation::cornerPoint(adjacency_, e, t);
    if (w >= 0) {
        for (const int c : adjacency_.pointCorners(w)) {
            if (adjacency_.cornerFace(c) == f) {
                u = static_cast<double>(c - c0);
                return true;
            }
        }
        return false;
    }
    for (const int c : adjacency_.edgeCorners(e)) {
        if (adjacency_.cornerFace(c) == f) {
            const float along = topology_.vertex(c) != adjacency_.edgeVertices(e)[0] ? 1.0f - t : t;
            const double value = static_cast<double>(c - c0) + static_cast<double>(along);
            u = value == static_cast<double>(k) ? 0.0 : value;
            return true;
        }
    }
    return false;
}

// Phase B.
void Cutter::countEntries()
{
    inside_ = BitArray(static_cast<std::size_t>(faceCount_));
    std::vector<int> z(static_cast<std::size_t>(faceCount_) + 1, 0);
    for (std::size_t s = 0; s < pieces_.size(); ++s) {
        const CurveSegment& piece = pieces_[s];
        const auto nF = static_cast<int>(piece.faces.size());
        if (nF == 0 || piece.size() < 2) {
            continue;
        }
        const bool closed = closed_[s];
        const int n = piece.size();
        for (int i = 0; i < nF && i + 1 < n; ++i) {
            const int f = piece.faces[static_cast<std::size_t>(i)];
            if (f < 0 || !needsCut(piece, i)) {
                continue;
            }
            const int e = piece.edges[static_cast<std::size_t>(i)];
            if (e < 0) {
                inside_.set(static_cast<std::size_t>(f));
            }
            const int before = i > 0 ? piece.faces[static_cast<std::size_t>(i) - 1] : (closed ? piece.faces.back() : -1);
            if ((i > 0 || closed) && before == f && e >= 0) {
                const int h = adjacency_.otherFace(e, f);
                if (h >= 0) {
                    ++z[static_cast<std::size_t>(h) + 1];
                }
            }
            ++z[static_cast<std::size_t>(f) + 1];
        }
        if (piece.edges.back() < 0 && piece.faces.back() >= 0) {
            inside_.set(static_cast<std::size_t>(piece.faces.back()));
        }
        if (!closed) {
            const std::vector<int>& kp = meshPoints_[s];
            for (const auto& [e, f, v] : {std::tuple{piece.edges.front(), piece.faces.front(), kp.front()},
                                          std::tuple{piece.edges.back(), piece.faces.back(), kp.back()}}) {
                if (v >= pointCount_ && e >= 0 && f >= 0) {
                    const int h = adjacency_.otherFace(e, f);
                    if (h >= 0) {
                        ++z[static_cast<std::size_t>(h) + 1];
                    }
                }
            }
        }
    }
    for (std::size_t f = 1; f < z.size(); ++f) {
        z[f] += z[f - 1];
    }
    entryStarts_ = std::move(z);
    entries_.assign(static_cast<std::size_t>(entryStarts_.back()), Entry{});
}

// Phase C.
void Cutter::writeEntries()
{
    std::vector<int> next(entryStarts_.begin(), entryStarts_.end() - 1);
    auto place = [&](int f) -> Entry& { return entries_[static_cast<std::size_t>(next[static_cast<std::size_t>(f)]++)]; };
    auto failed = [&](int f, Entry& entry) {
        const int c0 = topology_.vertex(topology_.faceBegin(f));
        entry = {0.0, 0.0, c0, c0};
    };
    auto pointEntry = [&](int h, float t, int e, int v) {
        Entry& entry = place(h);
        double u = 0.0;
        if (!perimeter(t, e, h, u)) {
            u = 0.0;
        }
        entry = {u, u, v, v};
    };
    for (std::size_t s = 0; s < pieces_.size(); ++s) {
        const CurveSegment& piece = pieces_[s];
        const auto nF = static_cast<int>(piece.faces.size());
        if (nF == 0 || piece.size() < 2) {
            continue;
        }
        const bool closed = closed_[s];
        const int n = piece.size();
        const std::vector<int>& kp = meshPoints_[s];
        for (int i = 0; i < nF && i + 1 < n; ++i) {
            const int f = piece.faces[static_cast<std::size_t>(i)];
            if (f < 0 || !needsCut(piece, i)) {
                continue;
            }
            const int e0 = piece.edges[static_cast<std::size_t>(i)];
            const int e1 = piece.edges[static_cast<std::size_t>(i) + 1];
            const float t0 = piece.t[static_cast<std::size_t>(i)];
            const float t1 = piece.t[static_cast<std::size_t>(i) + 1];
            const int before = i > 0 ? piece.faces[static_cast<std::size_t>(i) - 1] : (closed ? piece.faces.back() : -1);
            if ((i > 0 || closed) && before == f && e0 >= 0) {
                const int h = adjacency_.otherFace(e0, f);
                if (h >= 0) {
                    pointEntry(h, t0, e0, kp[static_cast<std::size_t>(i)]);
                }
            }
            Entry& entry = place(f);
            double tA = -1.0;
            double tB = -1.0;
            if (e0 >= 0 && !perimeter(t0, e0, f, tA)) {
                failed(f, entry);
                continue;
            }
            if (e1 >= 0 && !perimeter(t1, e1, f, tB)) {
                failed(f, entry);
                continue;
            }
            entry = {tA, tB, kp[static_cast<std::size_t>(i)], kp[static_cast<std::size_t>(i) + 1]};
            if (e1 >= 0 && tB < tA) {
                std::swap(entry.tA, entry.tB);
                std::swap(entry.vA, entry.vB);
            }
        }
        if (!closed) {
            const auto ends = {std::tuple{piece.edges.front(), piece.faces.front(), kp.front(), piece.t.front()},
                               std::tuple{piece.edges.back(), piece.faces.back(), kp.back(), piece.t.back()}};
            for (const auto& [e, f, v, t] : ends) {
                if (e >= 0 && f >= 0 && v >= pointCount_) {
                    const int h = adjacency_.otherFace(e, f);
                    if (h >= 0) {
                        pointEntry(h, t, e, v);
                    }
                }
            }
        }
    }
}

// Phase D.
void Cutter::mergePoints()
{
    const auto added = static_cast<int>(newPoints_.size());
    std::vector<int> rep(static_cast<std::size_t>(added));
    for (int i = 0; i < added; ++i) {
        rep[static_cast<std::size_t>(i)] = pointCount_ + i;
    }
    auto r = [&](int v) -> int& { return rep[static_cast<std::size_t>(v - pointCount_)]; };
    bool changed = false;
    std::vector<double> params;
    std::vector<int> border;
    std::vector<int> inner;
    std::vector<int> order;
    for (int f = 0; f < faceCount_; ++f) {
        const int begin = entryStarts_[static_cast<std::size_t>(f)];
        const int end = entryStarts_[static_cast<std::size_t>(f) + 1];
        if (begin == end) {
            continue;
        }
        params.clear();
        border.clear();
        inner.clear();
        for (int k = begin; k < end; ++k) {
            const Entry& entry = entries_[static_cast<std::size_t>(k)];
            for (const auto& [t, v] : {std::pair{entry.tA, entry.vA}, std::pair{entry.tB, entry.vB}}) {
                if (t >= 0.0) {
                    params.push_back(t);
                    border.push_back(v);
                } else if (std::find(inner.begin(), inner.end(), v) == inner.end()) {
                    inner.push_back(v);
                }
            }
        }
        order.resize(params.size());
        for (std::size_t i = 0; i < order.size(); ++i) {
            order[i] = static_cast<int>(i);
        }
        util::Introsort::sort(order.begin(), order.end(), [&](int x, int y) {
            return params[static_cast<std::size_t>(x)] < params[static_cast<std::size_t>(y)];
        });
        for (std::size_t j = 0; j + 1 < order.size(); ++j) {
            const auto x = static_cast<std::size_t>(order[j]);
            const auto y = static_cast<std::size_t>(order[j + 1]);
            const int bx = border[x];
            const int by = border[y];
            if (!(params[y] - params[x] > 1e-5) &&
                util::truncateToInt(params[x]) == util::truncateToInt(params[y]) && bx >= pointCount_ &&
                by >= pointCount_) {
                changed = true;
                if (r(bx) == bx) {
                    r(bx) = r(by);
                } else if (r(by) == by) {
                    r(by) = r(bx);
                }
            }
        }
        for (std::size_t a = 0; a + 1 < inner.size(); ++a) {
            for (std::size_t b = a + 1; b < inner.size(); ++b) {
                const int p = inner[a];
                const int q = inner[b];
                if (p < 0 || q < pointCount_ || r(q) != q) {
                    continue;
                }
                const Vec3f d = position(p) - position(q);
                const float d2 = dot(d, d);
                if (static_cast<double>(d2) < 1e-16) {
                    r(q) = r(p);
                    changed = true;
                }
            }
        }
    }
    if (!changed) {
        return;
    }
    for (int f = 0; f < faceCount_; ++f) {
        const int begin = entryStarts_[static_cast<std::size_t>(f)];
        const int end = entryStarts_[static_cast<std::size_t>(f) + 1];
        for (int k = begin; k < end; ++k) {
            Entry& entry = entries_[static_cast<std::size_t>(k)];
            for (int side = 0; side < 2; ++side) {
                double& t = side == 0 ? entry.tA : entry.tB;
                int& v = side == 0 ? entry.vA : entry.vB;
                if (v < pointCount_ || r(v) == v) {
                    continue;
                }
                const int target = r(v);
                if (t == -1.0) {
                    v = target;
                    continue;
                }
                for (int q = begin; q < end; ++q) {
                    const Entry& other = entries_[static_cast<std::size_t>(q)];
                    if (other.vA == target || other.vB == target) {
                        const double found = other.vA == target ? other.tA : other.tB;
                        if (found >= 0.0) {
                            t = found;
                            v = target;
                        }
                        break;
                    }
                }
            }
        }
    }
    for (auto& chain : chains_) {
        for (int& v : chain) {
            if (v >= pointCount_) {
                v = r(v);
            }
        }
    }
}

// Cleaning of the entries of a face without inside point (k corners). An
// entry from corner 0 to the closing side (k-1 -> 0) lies on that side.
void cleanEntries(std::vector<Entry>& entries, int k)
{
    std::size_t n = entries.size();
    for (std::size_t i = 0; i < n;) {
        Entry& e = entries[i];
        const double a = e.tA;
        const double b = e.tB;
        if (a == b) {
            if (e.vA != e.vB) {
                entries.erase(entries.begin() + static_cast<std::ptrdiff_t>(i));
                --n;
                continue;
            }
            ++i;
            continue;
        }
        const int ia = util::truncateToInt(a);
        const bool closingSide = a == 0.0 && util::truncateToInt(b) == k - 1;
        if (b > static_cast<double>(ia) + 1.0 && !closingSide) {
            ++i;
            continue;
        }
        const bool wa = isWhole(a);
        const bool wb = isWhole(b);
        if (wa && wb) {
            entries.erase(entries.begin() + static_cast<std::ptrdiff_t>(i));
            --n;
            continue;
        }
        if (wa) {
            e = {b, b, e.vB, e.vB};
        } else if (wb) {
            e = {a, a, e.vA, e.vA};
        } else {
            const Entry second{b, b, e.vB, e.vB};
            e = {a, a, e.vA, e.vA};
            entries.push_back(second);
        }
        ++i;
    }
}

// Splitting of one polygon along the chords of its entries (recursive).
class ChordSplitter {
public:
    ChordSplitter(Cutter& cutter, bool triangles) : cutter_(cutter), triangles_(triangles) {}

    void split(const Polygon& polygon, std::vector<Entry> entries, std::vector<Polygon>& out)
    {
        const int k = static_cast<int>(polygon.size());
        if (entries.empty()) {
            out.push_back(polygon);
            return;
        }
        Entry e0 = entries.front();
        normalize(e0, k);
        double a = 0.0;
        double b = 0.0;
        int i = 0;
        int j = 0;
        double fa = 0.0;
        double fb = 0.0;
        bool pointMode = false;
        Polygon p1;
        Polygon p2;
        int posB = 0;
        int posB2 = 0;
        while (true) {
            a = e0.tA;
            b = e0.tB;
            i = util::truncateToInt(a);
            j = util::truncateToInt(b);
            fa = a - i;
            fb = b - j;
            posB = fa != 0.0 ? i + 2 : i + 1;
            posB2 = fb != 0.0 ? j - i + 1 : j - i;
            const bool samePosition = a == b;
            if (fb == 0.0 && fa == 0.0 && j == i + 1) {
                // The chord is a side.
                restOrKeep(polygon, entries, out);
                return;
            }
            if (samePosition) {
                if (fa == 0.0) {
                    restOrKeep(polygon, entries, out);
                    return;
                }
                if (triangles_) {
                    // A point on a side of a triangle: chord to the opposite corner.
                    const int opposite = (i + 2) % k;
                    e0.vB = polygon[static_cast<std::size_t>(opposite)];
                    e0.tB = static_cast<double>(opposite);
                    if (0.0 <= e0.tB && e0.tB < e0.tA) {
                        std::swap(e0.tA, e0.tB);
                        std::swap(e0.vA, e0.vB);
                    }
                    continue;
                }
                pointMode = true;
                p1.assign(polygon.begin(), polygon.begin() + i + 1);
                p1.push_back(e0.vA);
                p1.insert(p1.end(), polygon.begin() + i + 1, polygon.end());
                posB = i + 1;
                break;
            }
            if (i == j) {
                // Both ends on one side: only this entry drops out.
                restOrKeep(polygon, entries, out);
                return;
            }
            p1.assign(polygon.begin(), polygon.begin() + i + 1);
            if (fa != 0.0) {
                p1.push_back(e0.vA);
            }
            p1.push_back(e0.vB);
            p1.insert(p1.end(), polygon.begin() + j + 1, polygon.end());
            p2.assign(1, e0.vA);
            p2.insert(p2.end(), polygon.begin() + i + 1, polygon.begin() + j + 1);
            if (fb != 0.0) {
                p2.push_back(e0.vB);
            }
            break;
        }

        std::vector<Entry> l1;
        std::vector<Entry> l2;
        std::vector<std::pair<double, int>> crossings;
        const auto region = [&](double z) { return z < a ? 1 : (z < b ? 2 : 3); };
        const int n1 = static_cast<int>(p1.size());
        const int n2 = static_cast<int>(p2.size());
        // Position of a border position z in its part polygon.
        const auto convert = [&](double z, int zone) {
            const int iz = util::truncateToInt(z);
            const double fz = z - iz;
            if (zone == 1) {
                return iz == i && fa != 0.0 ? static_cast<double>(i) + fz / fa : z;
            }
            if (zone == 2) {
                if (iz == i) {
                    return (fz - fa) / (1.0 - fa);
                }
                if (iz == j && fb != 0.0) {
                    return static_cast<double>(j - i) + fz / fb;
                }
                return z - i;
            }
            if (iz == j) {
                return static_cast<double>(posB) + (fz - fb) / (1.0 - fb);
            }
            return z + posB - j;
        };
        const auto addPart = [](std::vector<Entry>& list, Entry part, int corners) {
            if (part.tA == static_cast<double>(corners)) {
                part.tA = 0.0;
            }
            if (part.tB == static_cast<double>(corners)) {
                part.tB = 0.0;
            }
            if (part.tA != part.tB) {
                list.push_back(part);
            }
        };
        for (std::size_t m = 1; m < entries.size(); ++m) {
            Entry e = entries[m];
            normalize(e, k);
            if (e.vA == e0.vA && e.vB == e0.vB) {
                continue;
            }
            const double x = e.tA;
            const double y = e.tB;
            bool toSecond = false;
            double xp = 0.0;
            double yp = 0.0;
            bool done = false;
            // In point mode an end at the point entry lies on the inserted point.
            const auto atA = [&](bool second) { return second ? 0.0 : static_cast<double>(pointMode ? posB : posB - 1); };
            const auto atB = [&](bool second) { return second ? static_cast<double>(n2 - 1) : static_cast<double>(posB); };
            if (x == a || x == b) {
                toSecond = region(y) == 2;
                xp = x == a ? atA(toSecond) : atB(toSecond);
                yp = convert(y, region(y));
                done = true;
            } else if (y == a || y == b) {
                toSecond = region(x) == 2;
                yp = y == a ? atA(toSecond) : atB(toSecond);
                xp = convert(x, region(x));
                done = true;
            }
            if (!done) {
                const int rx = region(x);
                const int ry = region(y);
                const bool crossing = (rx == 2) != (ry == 2);
                xp = convert(x, rx);
                yp = convert(y, ry);
                if (crossing && !pointMode) {
                    const auto [s, w] = crossingPoint(e0, e, crossings);
                    if (rx == 2) {
                        addPart(l2, {xp, static_cast<double>(posB2) + (1.0 - s), e.vA, w}, n2);
                        addPart(l1, {static_cast<double>(posB - 1) + s, yp, w, e.vB}, n1);
                    } else {
                        addPart(l1, {xp, static_cast<double>(posB - 1) + s, e.vA, w}, n1);
                        addPart(l2, {static_cast<double>(posB2) + (1.0 - s), yp, w, e.vB}, n2);
                    }
                    continue;
                }
                toSecond = rx == 2 && ry == 2;
            }
            (toSecond ? l2 : l1).push_back({xp, yp, e.vA, e.vB});
        }
        if (pointMode) {
            if (n1 >= 3) {
                if (l1.empty()) {
                    out.push_back(p1);
                } else {
                    split(p1, std::move(l1), out);
                }
            }
            return;
        }
        if (n1 >= 3) {
            if (l1.empty()) {
                out.push_back(p1);
            } else {
                split(p1, std::move(l1), out);
            }
        }
        if (n2 >= 3) {
            if (l2.empty()) {
                out.push_back(p2);
            } else {
                split(p2, std::move(l2), out);
            }
        }
    }

private:
    static void normalize(Entry& e, int k)
    {
        if (e.tA == static_cast<double>(k)) {
            e.tA = 0.0;
        }
        if (e.tB == static_cast<double>(k)) {
            e.tB = 0.0;
        }
        if (0.0 <= e.tB && e.tB < e.tA) {
            std::swap(e.tA, e.tB);
            std::swap(e.vA, e.vB);
        }
    }

    // The first entry is done: split by the rest, or keep the polygon.
    void restOrKeep(const Polygon& polygon, std::vector<Entry>& entries, std::vector<Polygon>& out)
    {
        entries.erase(entries.begin());
        if (entries.empty()) {
            out.push_back(polygon);
        } else {
            split(polygon, std::move(entries), out);
        }
    }

    // Crossing of entry e with the chord e0: its parameter along e0 and its point.
    std::pair<double, int> crossingPoint(const Entry& e0, const Entry& e, std::vector<std::pair<double, int>>& known)
    {
        const Vec3d p = wide(cutter_.position(e0.vA));
        const Vec3d q = wide(cutter_.position(e0.vB));
        const Vec3d r = wide(cutter_.position(e.vA));
        const Vec3d sPoint = wide(cutter_.position(e.vB));
        Vec3d d1 = q - p;
        Vec3d d2 = sPoint - r;
        const double l1 = std::sqrt(dot(d1, d1));
        const double l2 = std::sqrt(dot(d2, d2));
        d1 = l1 < 1e-18 ? Vec3d{} : Vec3d{d1.x / l1, d1.y / l1, d1.z / l1};
        d2 = l2 < 1e-18 ? Vec3d{} : Vec3d{d2.x / l2, d2.y / l2, d2.z / l2};
        const double c = dot(d1, d2);
        const double den = 1.0 - c * c;
        int result = 3;
        double t1 = 0.0;
        Vec3d at{};
        if (!(den < 1e-50)) {
            const Vec3d w = r - p;
            t1 = (dot(d1, w) - c * dot(d2, w)) / den;
            const double t2 = (c * dot(d1, w) - dot(d2, w)) / den;
            at = {p.x + t1 * d1.x, p.y + t1 * d1.y, p.z + t1 * d1.z};
            result = (t1 < 0.0 || t1 > l1 || t2 < 0.0 || t2 > l2) ? 2 : 1;
        }
        double s = 0.0;
        const auto blend = [](const Vec3d& a, double wa, const Vec3d& b, double wb) {
            return Vec3d{a.x * wa + b.x * wb, a.y * wa + b.y * wb, a.z * wa + b.z * wb};
        };
        if (l1 == 0.0 || result != 1) {
            s = 0.5;
            at = {(p.x + q.x) * 0.5, (p.y + q.y) * 0.5, (p.z + q.z) * 0.5};
        } else {
            s = t1 / l1;
            const bool sharesA = e0.vA == e.vA || e0.vA == e.vB;
            const bool sharesB = e0.vB == e.vA || e0.vB == e.vB;
            if ((s == 0.0 && !sharesA) || s < 0.0) {
                s = 1e-12;
                at = blend(p, 0.999999999999, q, 1e-12);
            } else if ((s == 1.0 && !sharesB) || s > 1.0) {
                s = 0.999999999999;
                at = blend(p, 9.999778782798785e-13, q, 0.999999999999);
            }
        }
        for (const auto& [value, point] : known) {
            if (value == s) {
                return {s, point};
            }
        }
        const int point = cutter_.addPoint({static_cast<float>(at.x), static_cast<float>(at.y), static_cast<float>(at.z)});
        known.emplace_back(s, point);
        return {s, point};
    }

    Cutter& cutter_;
    bool triangles_;
};

// Moves a 2D point into the triangle (0, 1, 2) in up to three rounds.
void clampIntoTriangle(const std::array<CDT::V2d<double>, 3>& corners, CDT::V2d<double>& p)
{
    for (int round = 0; round < 3; ++round) {
        bool outside = false;
        for (int side = 0; side < 3; ++side) {
            const CDT::V2d<double>& a = corners[static_cast<std::size_t>(side)];
            const CDT::V2d<double>& b = corners[static_cast<std::size_t>((side + 1) % 3)];
            const CDT::V2d<double>& o = corners[static_cast<std::size_t>((side + 2) % 3)];
            const double ux = b.x - a.x;
            const double uy = b.y - a.y;
            const double crossed = ux * (p.y - a.y) - uy * (p.x - a.x);
            if (crossed < 0.0) {
                outside = true;
                const double len2 = ux * ux + uy * uy;
                double t = len2 > 0.0 ? ((p.x - a.x) * ux + (p.y - a.y) * uy) / len2 : 0.0;
                t = std::clamp(t, 0.0, 1.0);
                const double fx = a.x + t * ux;
                const double fy = a.y + t * uy;
                p.x = fx + 0.001 * (o.x - fx);
                p.y = fy + 0.001 * (o.y - fy);
            }
        }
        if (!outside) {
            return;
        }
    }
}

// Constrained Delaunay triangulation of a triangle with its entries.
bool triangulateWithEntries(Cutter& cutter, Polygon& polygon, std::span<const Entry> entries)
{
    if (polygon.size() != 3) {
        return false; // only triangles carry inside points
    }
    const Vec3f o = cutter.position(polygon[0]);
    const Vec3f d1 = cutter.position(polygon[1]) - o;
    const Vec3f d2 = cutter.position(polygon[2]) - o;
    Vec3d x = wide(d1);
    const double len = normalizeInPlace(x);
    const Vec3d e2 = wide(d2);
    Vec3d y = e2 - Vec3d{x.x * dot(e2, x), x.y * dot(e2, x), x.z * dot(e2, x)};
    normalizeInPlace(y);
    std::vector<CDT::V2d<double>> points;
    std::vector<int> meshPoint;
    const std::array<CDT::V2d<double>, 3> corners{CDT::V2d<double>::make(0.0, 0.0), CDT::V2d<double>::make(len, 0.0),
                                                  CDT::V2d<double>::make(dot(x, e2), dot(y, e2))};
    for (int i = 0; i < 3; ++i) {
        points.push_back(corners[static_cast<std::size_t>(i)]);
        meshPoint.push_back(polygon[static_cast<std::size_t>(i)]);
    }
    std::vector<std::pair<double, int>> border{{0.0, 0}, {1.0, 1}, {2.0, 2}};
    std::vector<CDT::Edge> constraints;
    auto indexOf = [&](int v, double t) {
        for (std::size_t i = 0; i < meshPoint.size(); ++i) {
            if (meshPoint[i] == v) {
                return static_cast<int>(i);
            }
        }
        const Vec3d d = wide(cutter.position(v) - o);
        CDT::V2d<double> p = CDT::V2d<double>::make(dot(x, d), dot(y, d));
        if (t == -1.0) {
            clampIntoTriangle(corners, p);
        }
        points.push_back(p);
        meshPoint.push_back(v);
        return static_cast<int>(points.size()) - 1;
    };
    for (const Entry& entry : entries) {
        const int a = indexOf(entry.vA, entry.tA);
        const int b = indexOf(entry.vB, entry.tB);
        if (a != b) {
            constraints.emplace_back(static_cast<CDT::VertInd>(a), static_cast<CDT::VertInd>(b));
        }
        if (entry.tA >= 0.0) {
            border.emplace_back(entry.tA, a);
        }
        if (entry.tB >= 0.0) {
            border.emplace_back(entry.tB, b);
        }
    }
    util::Introsort::sort(border.begin(), border.end(), [](const auto& l, const auto& r) { return l.first < r.first; });
    for (std::size_t i = 0; i < border.size(); ++i) {
        const auto& previous = border[(i + border.size() - 1) % border.size()];
        const auto& current = border[i];
        if (previous.first != current.first && previous.second != current.second) {
            constraints.emplace_back(static_cast<CDT::VertInd>(previous.second), static_cast<CDT::VertInd>(current.second));
        }
    }
    // Equal 2D points stop the triangulation; the face stays as it is.
    std::map<std::pair<double, double>, int> seen;
    for (const auto& p : points) {
        if (!seen.emplace(std::pair{p.x, p.y}, 0).second) {
            return false;
        }
    }
    try {
        CDT::Triangulation<double> cdt(CDT::VertexInsertionOrder::Auto, CDT::IntersectingConstraintEdges::Resolve, 0.0);
        cdt.insertVertices(points);
        cdt.insertEdges(constraints);
        cdt.eraseOuterTriangles();
        for (std::size_t i = points.size(); i < cdt.vertices.size(); ++i) {
            const auto& v = cdt.vertices[i];
            const Vec3f p{static_cast<float>(v.x * x.x + v.y * y.x) + o.x, static_cast<float>(v.x * x.y + v.y * y.y) + o.y,
                          static_cast<float>(v.x * x.z + v.y * y.z) + o.z};
            meshPoint.push_back(cutter.addPoint(p));
        }
        Polygon triangles;
        for (const CDT::Triangle& t : cdt.triangles) {
            for (const auto v : t.vertices) {
                triangles.push_back(meshPoint[static_cast<std::size_t>(v)]);
            }
        }
        polygon = std::move(triangles);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

// Phase E.
void Cutter::rebuildFaces()
{
    removed_ = BitArray(static_cast<std::size_t>(faceCount_));
    const bool triangles = topology_.mode() == FaceMode::Triangles;
    ChordSplitter splitter(*this, triangles);
    std::vector<Polygon> parts;
    for (int f = 0; f < faceCount_; ++f) {
        const int begin = entryStarts_[static_cast<std::size_t>(f)];
        const int end = entryStarts_[static_cast<std::size_t>(f) + 1];
        if (begin == end) {
            continue;
        }
        Polygon polygon;
        for (int c = topology_.faceBegin(f); c < topology_.faceEnd(f); ++c) {
            polygon.push_back(topology_.vertex(c));
        }
        parts.clear();
        const std::span<const Entry> own(entries_.data() + begin, static_cast<std::size_t>(end - begin));
        if (inside_.test(static_cast<std::size_t>(f))) {
            Polygon result = polygon;
            if (triangulateWithEntries(*this, result, own)) {
                for (std::size_t k = 0; k + 2 < result.size(); k += 3) {
                    parts.push_back({result[k], result[k + 1], result[k + 2]});
                }
            } else {
                parts.push_back(polygon);
            }
        } else {
            std::vector<Entry> list(own.begin(), own.end());
            cleanEntries(list, static_cast<int>(polygon.size()));
            if (list.empty()) {
                continue;
            }
            splitter.split(polygon, std::move(list), parts);
        }
        removed_.set(static_cast<std::size_t>(f));
        for (Polygon& part : parts) {
            polygons_.push_back(std::move(part));
            origins_.push_back(f);
        }
        // Every new point of this face gets its own source row.
        while (sources_.rows() < static_cast<int>(newPoints_.size())) {
            addFaceRow(f);
        }
    }
}

// Phase F: the new polygons of a triangle mesh as triangles; other meshes
// take them as they are.
Mesh Cutter::triangulatedFaces(std::vector<int>& origin) const
{
    Mesh added;
    std::vector<int> corners;
    std::vector<int> starts{0};
    origin.clear();
    const bool triangles = topology_.mode() == FaceMode::Triangles;
    std::vector<Vec3f> polygon;
    for (std::size_t p = 0; p < polygons_.size(); ++p) {
        const Polygon& face = polygons_[p];
        const auto k = face.size();
        auto emit = [&](int a, int b, int c) {
            corners.insert(corners.end(), {a, b, c});
            starts.push_back(static_cast<int>(corners.size()));
            origin.push_back(origins_[p]);
        };
        if (!triangles || k == 3) {
            corners.insert(corners.end(), face.begin(), face.end());
            starts.push_back(static_cast<int>(corners.size()));
            origin.push_back(origins_[p]);
        } else if (k == 4) {
            const Vec3f d02 = position(face[0]) - position(face[2]);
            const Vec3f d13 = position(face[1]) - position(face[3]);
            if (dot(d02, d02) <= dot(d13, d13)) {
                emit(face[0], face[1], face[2]);
                emit(face[0], face[2], face[3]);
            } else {
                emit(face[1], face[2], face[3]);
                emit(face[1], face[3], face[0]);
            }
        } else if (k > 4) {
            polygon.clear();
            for (const int v : face) {
                polygon.push_back(position(v));
            }
            for (const auto& t : prep::PolygonTriangulator::triangulate(polygon)) {
                emit(face[static_cast<std::size_t>(t[0])], face[static_cast<std::size_t>(t[1])],
                     face[static_cast<std::size_t>(t[2])]);
            }
        }
    }
    std::vector<int> cornerList = std::move(corners);
    added.editTopology().swapPolygons(starts, cornerList);
    if (triangles) {
        added.editTopology().trySetUniformMode(FaceMode::Triangles);
    } else if (topology_.mode() == FaceMode::Quads) {
        added.editTopology().trySetUniformMode(FaceMode::Quads);
    }
    added.topologyChanged();
    return added;
}

// Phase G: only a real change is applied.
bool Cutter::apply()
{
    if (newPoints_.empty() && !removed_.any()) {
        return false;
    }
    std::vector<int> origin;
    Mesh added = triangulatedFaces(origin);
    Modification change;
    change.appendPoints = newPoints_;
    change.pointSources = sources_;
    change.pointWeights = weights_;
    change.deleteFaces = removed_;
    change.appendFaces = &added;
    change.points.state = 2; // old points keep their numbers
    change.faces.newToOld = std::move(origin);
    return MeshModifier::apply(mesh_, change);
}

// Phase H.
void Cutter::repairChains()
{
    const Adjacency& adjacency = mesh_.adjacency();
    for (std::size_t s = 0; s < chains_.size(); ++s) {
        std::vector<int>& chain = chains_[s];
        std::vector<int>& numbers = numbers_[s];
        for (std::size_t j = 1; j < chain.size();) {
            const int a = chain[j - 1];
            const int b = chain[j];
            if (adjacency.edgeBetween(a, b) >= 0) {
                ++j;
                continue;
            }
            int best = -1;
            float bestLength = FLT_MAX;
            const auto neighbours = adjacency.pointNeighbours(a);
            for (const int w : neighbours) {
                if (w >= pointCount_ && adjacency.edgeBetween(w, b) >= 0) {
                    const float d = distance(mesh_.points()[static_cast<std::size_t>(b)], mesh_.points()[static_cast<std::size_t>(w)]) +
                                    distance(mesh_.points()[static_cast<std::size_t>(w)], mesh_.points()[static_cast<std::size_t>(a)]);
                    if (d < bestLength) {
                        bestLength = d;
                        best = w;
                    }
                }
            }
            const int number = numbers[j - 1];
            if (best >= 0) {
                chain.insert(chain.begin() + static_cast<std::ptrdiff_t>(j), best);
                numbers.insert(numbers.begin() + static_cast<std::ptrdiff_t>(j), number);
                j += 2;
                continue;
            }
            // The shortest edge path brings the points between a and b.
            std::vector<int> path = CurveProjector::edgePath(mesh_, a, b);
            if (path.size() >= 2) {
                path.pop_back();
                chain.insert(chain.begin() + static_cast<std::ptrdiff_t>(j), path.begin(), path.end());
                numbers.insert(numbers.begin() + static_cast<std::ptrdiff_t>(j), path.size(), number);
                j += path.size();
            }
            ++j;
        }
    }
}

SegmentCutter::Result Cutter::run()
{
    SegmentCutter::Result result;
    if (pieces_.empty()) {
        return result;
    }
    toMeshPoints();
    countEntries();
    writeEntries();
    mergePoints();
    rebuildFaces();
    result.changed = apply();
    repairChains();
    if (options_.compactPoints) {
        std::vector<int> renumber;
        if (prep::MeshRepair::removeUnusedPoints(mesh_, &renumber) > 0) {
            result.changed = true;
            for (auto& chain : chains_) {
                for (int& v : chain) {
                    v = v >= 0 ? renumber[static_cast<std::size_t>(v)] : v;
                }
            }
        }
    }
    result.chains = std::move(chains_);
    result.numbers = std::move(numbers_);
    return result;
}

} // namespace

SegmentCutter::Result SegmentCutter::cut(Mesh& mesh, std::span<const CurveSegment> pieces, const Options& options)
{
    return Cutter(mesh, pieces, options).run();
}

} // namespace areko::remesher::curves
