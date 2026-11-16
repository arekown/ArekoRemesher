#include "segment_preparation.h"

#include <algorithm>
#include <cstddef>

namespace areko::remesher::curves {

using namespace areko::remesher::mesh;

namespace {

constexpr int pullRounds = 5;

bool atCorner(float t) { return t == 0.0f || t == 1.0f; }

void erasePoint(CurveSegment& piece, int i)
{
    const auto k = static_cast<std::ptrdiff_t>(i);
    piece.points.erase(piece.points.begin() + k);
    piece.edges.erase(piece.edges.begin() + k);
    piece.t.erase(piece.t.begin() + k);
    if (static_cast<std::size_t>(i) < piece.sources.size()) {
        piece.sources.erase(piece.sources.begin() + k);
    }
}

} // namespace

int SegmentPreparation::cornerPoint(const Adjacency& adjacency, int edge, float t)
{
    if (edge < 0) {
        return -1;
    }
    if (t == 0.0f) {
        return adjacency.edgeVertices(edge)[0];
    }
    if (t == 1.0f) {
        return adjacency.edgeVertices(edge)[1];
    }
    return -1;
}

bool SegmentPreparation::isClosed(const Adjacency& adjacency, const CurveSegment& piece)
{
    const int n = piece.size();
    if (n < 3) {
        return false;
    }
    const int e0 = piece.edges.front();
    const int e1 = piece.edges.back();
    if (e0 < 0 || e1 < 0) {
        return false;
    }
    const float t0 = piece.t.front();
    const float t1 = piece.t.back();
    if (!atCorner(t0)) {
        return t0 == t1 && e0 == e1;
    }
    return atCorner(t1) && cornerPoint(adjacency, e0, t0) == cornerPoint(adjacency, e1, t1);
}

std::vector<CurveSegment> SegmentPreparation::prepare(Mesh& mesh, std::span<const CurveSegment* const> pieces,
                                                      const Options& options)
{
    std::vector<CurveSegment> out;
    out.reserve(pieces.size());
    const Adjacency& adjacency = mesh.adjacency();
    std::vector<int> movedPoints;
    std::vector<Vec3f> movedTo;
    const float tau = options.snap;
    for (const CurveSegment* source : pieces) {
        CurveSegment piece = *source;
        // A: snap parameters close to an end point onto it.
        for (int i = 0; i < piece.size(); ++i) {
            const int e = piece.edges[static_cast<std::size_t>(i)];
            if (e < 0) {
                continue;
            }
            float& t = piece.t[static_cast<std::size_t>(i)];
            if (atCorner(t) || (tau <= t && t <= 1.0f - tau)) {
                continue;
            }
            const auto ends = adjacency.edgeVertices(e);
            int v = -1;
            if (t < tau) {
                t = 0.0f;
                v = ends[0];
            } else if (t > 1.0f - tau) {
                t = 1.0f;
                v = ends[1];
            }
            if (v >= 0) {
                movedPoints.push_back(v);
                movedTo.push_back(piece.points[static_cast<std::size_t>(i)]);
                if (options.snapToPoint) {
                    piece.points[static_cast<std::size_t>(i)] = mesh.points()[static_cast<std::size_t>(v)];
                }
            }
        }
        // B: points between two corner points on the edge of these points.
        int last = -1;
        for (int j = 0; j < piece.size(); ++j) {
            const int vj = cornerPoint(adjacency, piece.edges[static_cast<std::size_t>(j)], piece.t[static_cast<std::size_t>(j)]);
            if (vj < 0) {
                continue;
            }
            if (last >= 0 && last + 1 < j) {
                const int vi = cornerPoint(adjacency, piece.edges[static_cast<std::size_t>(last)],
                                           piece.t[static_cast<std::size_t>(last)]);
                const int e = adjacency.edgeBetween(vi, vj);
                bool along = e >= 0;
                for (int k = last + 1; along && k < j; ++k) {
                    along = piece.edges[static_cast<std::size_t>(k)] == e;
                }
                if (along) {
                    for (int k = last + 1; k < j; ++k) {
                        piece.t[static_cast<std::size_t>(k)] = 0.0f;
                    }
                }
            }
            last = j;
        }
        // C: zig-zags over one edge. A closed piece runs round: b and c wrap
        // to its start, and the face before entry 0 is its last face.
        const bool closed = isClosed(adjacency, piece);
        int n = piece.size();
        int m = closed ? n - 1 : n;
        int limit = closed ? n - 2 : n - 3;
        for (int i = closed ? 0 : 1; i <= limit;) {
            const auto faces = static_cast<int>(piece.faces.size());
            const int b = m > 0 ? (i + 1) % m : 0;
            const int c = m > 0 ? (i + 2) % m : 0;
            const int before = i > 0 ? piece.faces[static_cast<std::size_t>(i - 1)] : (closed && faces > 0 ? piece.faces.back() : -1);
            bool zigzag = before >= 0 && (closed || c != 0) && i < faces && b < faces && c < faces;
            if (zigzag) {
                const int e = piece.edges[static_cast<std::size_t>(i)];
                zigzag = e == piece.edges[static_cast<std::size_t>(b)] && e == piece.edges[static_cast<std::size_t>(c)] &&
                         before == piece.faces[static_cast<std::size_t>(b)] &&
                         piece.faces[static_cast<std::size_t>(i)] == piece.faces[static_cast<std::size_t>(c)];
            }
            if (!zigzag) {
                ++i;
                continue;
            }
            float& t = piece.t[static_cast<std::size_t>(i)];
            t = ((t + piece.t[static_cast<std::size_t>(b)]) + piece.t[static_cast<std::size_t>(c)]) / 3.0f;
            Vec3f& p = piece.points[static_cast<std::size_t>(i)];
            const Vec3f sum = (p + piece.points[static_cast<std::size_t>(b)]) + piece.points[static_cast<std::size_t>(c)];
            p = {sum.x / 3.0f, sum.y / 3.0f, sum.z / 3.0f};
            // Entries b and c and the faces of sections i and b go; the later
            // index first, so the earlier one stays valid.
            erasePoint(piece, std::max(b, c));
            erasePoint(piece, std::min(b, c));
            piece.faces.erase(piece.faces.begin() + std::max(i, b));
            piece.faces.erase(piece.faces.begin() + std::min(i, b));
            n -= 2;
            m -= 2;
            limit -= 2;
            i = std::max(i - 2, 0) + 1;
        }
        // D: neighbours of corner points that lie on an edge of the point.
        n = piece.size();
        for (int round = 0; round < pullRounds; ++round) {
            bool changed = false;
            const int from = closed ? 0 : 1;
            const int to = closed ? n - 1 : n - 2;
            for (int j = from; j <= to && n > 0; ++j) {
                const int vj = cornerPoint(adjacency, piece.edges[static_cast<std::size_t>(j)], piece.t[static_cast<std::size_t>(j)]);
                if (vj < 0) {
                    continue;
                }
                const int before = (j - 1 + n) % n;
                const int after = (j + 1) % n;
                // Both neighbours need an edge, else neither is pulled.
                if (piece.edges[static_cast<std::size_t>(before)] < 0 || piece.edges[static_cast<std::size_t>(after)] < 0) {
                    continue;
                }
                for (const int k : {before, after}) {
                    const int e = piece.edges[static_cast<std::size_t>(k)];
                    float& t = piece.t[static_cast<std::size_t>(k)];
                    if (atCorner(t)) {
                        continue;
                    }
                    const auto ends = adjacency.edgeVertices(e);
                    if (ends[0] == vj) {
                        t = 0.0f;
                        changed = true;
                    } else if (ends[1] == vj) {
                        t = 1.0f;
                        changed = true;
                    }
                }
            }
            if (!changed) {
                break;
            }
        }
        out.push_back(std::move(piece));
    }
    if (!movedPoints.empty()) {
        mesh.movePoints(movedPoints, movedTo);
    }
    return out;
}

} // namespace areko::remesher::curves
