#include "curve_cut_transfer.h"

#include "curves/curve_projector.h"
#include "curves/segment_cutter.h"
#include "curves/segment_preparation.h"

#include <utility>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;
using curves::CurveSegment;

namespace {

// Moves the entries between two reached polyline points of a piece onto the
// crease chain of the triangle mesh: by their share of the arc length between
// the two reached points they get a place on the chain, which gives their
// position and their source. Every such entry gets the code -2 - (number of
// its source in `sources`).
class ChainResampler {
public:
    ChainResampler(const Mesh& triangles, std::span<const int> chain, std::span<const float> keys,
                   std::vector<PointSource>& sources)
        : triangles_(triangles), chain_(chain), keys_(keys), sources_(sources)
    {
    }

    // Every pair of consecutive reached entries (code >= 0) of `piece`; the
    // entries before the first and after the last one stay as they are.
    void anchor(CurveSegment& piece)
    {
        int a = -1;
        for (int i = 0; i < static_cast<int>(piece.points.size()); ++i) {
            if (piece.sources[static_cast<std::size_t>(i)] >= 0) {
                if (a != -1) {
                    section(piece, a, i);
                }
                a = i;
            }
        }
    }

private:
    void section(CurveSegment& piece, int a, int i)
    {
        const auto count = static_cast<std::size_t>(i - a);
        share_.assign(count + 1, 0.0f);
        for (std::size_t k = 1; k <= count; ++k) {
            const auto at = static_cast<std::size_t>(a) + k;
            share_[k] = length(piece.points[at] - piece.points[at - 1]) + share_[k - 1];
        }
        const float total = share_[count];
        for (std::size_t r = 1; r <= count; ++r) {
            share_[r] = share_[r] / total;
        }
        const float sa = keys_[static_cast<std::size_t>(piece.sources[static_cast<std::size_t>(a)])];
        const float si = keys_[static_cast<std::size_t>(piece.sources[static_cast<std::size_t>(i)])];
        // A closed chain runs round when the second place lies before the first.
        const float last = static_cast<float>(static_cast<int>(chain_.size()) - 1);
        const bool round = sa > si;
        const float span = round ? (si + last) - sa : si - sa;
        for (int k = a + 1; k < i; ++k) {
            float u = span * share_[static_cast<std::size_t>(k - a)] + sa;
            if (round && u >= last) {
                u = u - last;
            }
            const auto at = static_cast<std::size_t>(k);
            piece.points[at] = chainPoint(u, piece.points[at]);
            piece.sources[at] = -2 - static_cast<int>(sources_.size());
            sources_.push_back(found_);
        }
    }

    // Position at place u of the chain; its source goes to found_. A place
    // inside chain edge j gives the edge (none when the two points have no
    // edge between them), a place at or behind the last point gives that
    // point. A place before the chain (not a number) keeps `position` and
    // gives no source.
    Vec3f chainPoint(float u, Vec3f position)
    {
        const auto& points = triangles_.points();
        const int n = static_cast<int>(chain_.size());
        const int j = truncToInt(u);
        if (j >= n - 1) {
            const int v = chain_[static_cast<std::size_t>(n - 1)];
            found_ = {1, v};
            return points[static_cast<std::size_t>(v)];
        }
        if (j < 0) {
            found_ = {0, -1};
            return position;
        }
        const int v0 = chain_[static_cast<std::size_t>(j)];
        const int v1 = chain_[static_cast<std::size_t>(j) + 1];
        const int e = triangles_.adjacency().edgeBetween(v0, v1);
        found_ = e >= 0 ? PointSource{2, e} : PointSource{0, -1};
        const float f = u - static_cast<float>(j);
        const float g = 1.0f - f;
        const Vec3f p0 = points[static_cast<std::size_t>(v0)];
        const Vec3f p1 = points[static_cast<std::size_t>(v1)];
        return {f * p1.x + g * p0.x, f * p1.y + g * p0.y, f * p1.z + g * p0.z};
    }

    const Mesh& triangles_;
    std::span<const int> chain_;
    std::span<const float> keys_;
    std::vector<PointSource>& sources_;
    std::vector<float> share_;
    PointSource found_;
};

} // namespace

void CurveCutTransfer::markEdges(const Adjacency& adjacency, std::span<const int> chain, BitArray& edges)
{
    for (std::size_t j = 0; j + 1 < chain.size(); ++j) {
        // A point twice in a row joins nothing.
        const int e = adjacency.edgeBetween(chain[j], chain[j + 1]);
        if (e >= 0) {
            edges.set(static_cast<std::size_t>(e));
        }
    }
}

void CurveCutTransfer::run(const Mesh& triangles, const CurveMarks& marks, std::span<const CurvePolyline> polylines,
                           RawMesh& raw)
{
    Mesh& mesh = raw.mesh;
    const curves::CurveProjector::Options projection{.tolerance = 5e-6f, .chordCheck = true, .sourceCodes = true};
    std::vector<curves::LocatedPolyline> located(polylines.size());
    for (std::size_t c = 0; c < polylines.size(); ++c) {
        for (const int v : polylines[c].points) {
            located[c].points.push_back(mesh.points()[static_cast<std::size_t>(v)]);
            located[c].locations.push_back({curves::Location::Point, v});
        }
    }
    std::vector<std::vector<CurveSegment>> projected = curves::CurveProjector::projectAll(mesh, located, projection);
    std::vector<PointSource> sources;
    std::vector<CurveSegment> pieces;
    for (std::size_t c = 0; c < polylines.size(); ++c) {
        ChainResampler resampler(triangles, marks.chainPoints.row(static_cast<int>(c)), polylines[c].keys, sources);
        for (CurveSegment& piece : projected[c]) {
            resampler.anchor(piece);
            pieces.push_back(std::move(piece));
        }
    }
    std::vector<const CurveSegment*> list;
    list.reserve(pieces.size());
    for (const CurveSegment& piece : pieces) {
        list.push_back(&piece);
    }
    const std::vector<CurveSegment> prepared =
        curves::SegmentPreparation::prepare(mesh, list, {.snap = 0.0f, .snapToPoint = false});
    const int before = mesh.pointCount();
    const curves::SegmentCutter::Result result =
        curves::SegmentCutter::cut(mesh, prepared, {.compactPoints = false, .moveCornerPoints = true});

    // The point data follow the new points; every new point lies on a curve.
    OutputPoints& points = raw.points;
    const auto count = static_cast<std::size_t>(mesh.pointCount());
    points.resize(count, PointSource{0, -1}, 0.0f);
    for (auto v = static_cast<std::size_t>(before); v < count; ++v) {
        points.onCurve.set(v);
    }
    for (std::size_t s = 0; s < result.chains.size(); ++s) {
        const std::vector<int>& chain = result.chains[s];
        const std::vector<int>& numbers = result.numbers[s];
        const std::vector<int>& codes = prepared[s].sources;
        for (std::size_t j = 0; j < chain.size() && j < numbers.size(); ++j) {
            const int v = chain[j];
            const int k = numbers[j];
            if (v >= before && k >= 0 && static_cast<std::size_t>(k) < codes.size()) {
                const int code = codes[static_cast<std::size_t>(k)];
                if (code < -1) {
                    points.source[static_cast<std::size_t>(v)] = sources[static_cast<std::size_t>(-2 - code)];
                }
            }
        }
    }
    const Adjacency& adjacency = mesh.adjacency();
    raw.curveEdges = BitArray(static_cast<std::size_t>(adjacency.edgeCount()));
    for (const std::vector<int>& chain : result.chains) {
        markEdges(adjacency, chain, raw.curveEdges);
    }
}

} // namespace areko::remesher::extract
