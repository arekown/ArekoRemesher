#include "hole_closing.h"

#include "mesh/boundary_loops.h"
#include "mesh/connected_parts.h"
#include "mesh/modification.h"
#include "prep/mesh_repair.h"

#include <algorithm>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

constexpr std::size_t largestHole = 1000; // loops from this size on stay open

void appendRow(Csr& csr, std::span<const int> row)
{
    csr.entries.insert(csr.entries.end(), row.begin(), row.end());
    csr.starts.push_back(static_cast<int>(csr.entries.size()));
}

} // namespace

HoleClosing::HoleClosing(Mesh& output, OutputPoints& points, const Mesh& original, int quadBudget,
                         const FrozenBorderInput* frozen, Options options)
    : output_(output), points_(points), original_(original), quadBudget_(quadBudget), frozen_(frozen), options_(options)
{
}

void HoleClosing::run()
{
    splitFans();
    const BoundaryLoops loops(output_);
    if (options_.rebuildPolygons || frozen_ != nullptr) {
        closeLoops(loops);
        if (frozen_ == nullptr) {
            if (options_.followBorders) {
                closeNotches(loops);
            }
        } else if (reachingLoops_.rows() > 0) {
            closeGapRings();
        }
    }
    apply();
}

// a) Border points outside every border loop and the end points of edges
// with more than two corners get one point per further fan.
void HoleClosing::splitFans()
{
    const int oldCount = output_.pointCount();
    BitArray mask(static_cast<std::size_t>(oldCount));
    {
        const Adjacency& adjacency = output_.adjacency();
        for (const int v : adjacency.boundaryPoints()) {
            mask.set(static_cast<std::size_t>(v));
        }
        const BoundaryLoops loops(output_);
        for (const int v : loops.points().entries) {
            mask.set(static_cast<std::size_t>(v), false);
        }
        for (const int e : adjacency.nonManifoldEdges()) {
            for (const int v : adjacency.edgeVertices(e)) {
                mask.set(static_cast<std::size_t>(v));
            }
        }
    }
    std::vector<int> sources;
    prep::MeshRepair::splitFans(output_, mask, sources);
    const int newCount = output_.pointCount();
    if (newCount == oldCount) {
        return;
    }
    points_.resize(static_cast<std::size_t>(newCount), PointSource{}, 1.0f);
    for (std::size_t j = 0; j < sources.size(); ++j) {
        const auto from = static_cast<std::size_t>(sources[j]);
        const std::size_t to = static_cast<std::size_t>(oldCount) + j;
        points_.onBorderEdge.set(to, points_.onBorderEdge.testBounded(from));
        points_.reachedBorder.set(to, points_.reachedBorder.testBounded(from));
        points_.source[to] = points_.source[from];
        points_.weight[to] = points_.weight[from];
    }
}

// b) Holes: loops that reach the border are kept for the gap rings, the
// others become polygons or close by a join.
void HoleClosing::closeLoops(const BoundaryLoops& loops)
{
    const float threshold = original_.triangleGeometry().totalArea() / (static_cast<float>(quadBudget_) * 500.0f);
    std::vector<int> hole;
    for (int r = 0; r < loops.count(); ++r) {
        const auto loop = loops.points().row(r);
        const bool reaches = std::any_of(loop.begin(), loop.end(), [&](int v) {
            return points_.reachedBorder.testBounded(static_cast<std::size_t>(v));
        });
        if (reaches) {
            if (frozen_ != nullptr) {
                appendRow(reachingLoops_, loop);
            }
            continue;
        }
        if (!options_.rebuildPolygons || loop.size() >= largestHole) {
            continue;
        }
        hole.assign(loop.rbegin(), loop.rend());
        if (hole.size() == 4 && collapse(hole, threshold)) {
            continue;
        }
        if (hole.size() >= 3) {
            addPolygon(hole);
        }
    }
}

// A thin four point hole: all corners of q0 go to q2 or of q1 to q3.
bool HoleClosing::collapse(std::span<const int> quad, float threshold)
{
    const auto& p = output_.points();
    auto at = [&](int k) { return p[static_cast<std::size_t>(quad[static_cast<std::size_t>(k)])]; };
    const Vec3f diagonal02 = at(0) - at(2);
    const Vec3f diagonal13 = at(1) - at(3);
    const float d02 = dot(diagonal02, diagonal02);
    const float d13 = dot(diagonal13, diagonal13);
    int from = -1;
    int to = -1;
    if (d02 < d13 / 9.0f && d02 < threshold) {
        from = quad[0];
        to = quad[2];
    } else if (d13 < d02 / 9.0f && d13 < threshold) {
        from = quad[1];
        to = quad[3];
    } else {
        return false;
    }
    for (const int c : output_.adjacency().pointCorners(from)) {
        joinedCorners_.push_back(c);
        joinTargets_.push_back(to);
    }
    return true;
}

// c) Without frozen borders: the stretch between two points on border edges
// closes against the loop direction.
void HoleClosing::closeNotches(const BoundaryLoops& loops)
{
    std::vector<int> run;
    for (int r = 0; r < loops.count(); ++r) {
        const auto loop = loops.points().row(r);
        const std::size_t n = loop.size();
        std::size_t s = 0;
        while (s < n && !points_.onBorderEdge.testBounded(static_cast<std::size_t>(loop[s]))) {
            ++s;
        }
        if (s == n) {
            continue;
        }
        run.assign(1, loop[s]);
        for (std::size_t i = s; i <= s + n; ++i) {
            const int q = loop[i % n];
            if (!points_.onBorderEdge.testBounded(static_cast<std::size_t>(q))) {
                run.push_back(q);
                continue;
            }
            if (run.size() >= 2 && run.size() != n) {
                run.push_back(q);
                std::reverse(run.begin(), run.end());
                if (run.size() >= 3) {
                    addPolygon(run);
                }
            }
            run.assign(1, q);
        }
    }
}

// c) With frozen borders: the rings between the loops that reach the border
// and the kept border chains.
void HoleClosing::closeGapRings()
{
    const BoundaryLoops originalLoops(original_);
    std::vector<int> lid;
    for (int r = 0; r < originalLoops.count(); ++r) {
        lid.clear();
        for (const int v : originalLoops.points().row(r)) {
            if (frozen_->marks[static_cast<std::size_t>(v)].value == 1) {
                lid.push_back(v);
            }
        }
        std::reverse(lid.begin(), lid.end());
        for (int& v : lid) {
            v = frozen_->originalToOutput[static_cast<std::size_t>(v)];
        }
        if (!lid.empty()) {
            appendRow(reachingLoops_, lid);
        }
    }

    // The rings are the border loops of the mesh of all these polygons.
    Topology rings;
    std::vector<int> starts = reachingLoops_.starts;
    std::vector<int> corners = reachingLoops_.entries;
    rings.swapPolygons(starts, corners);
    const Adjacency adjacency(rings, output_.pointCount());
    const ConnectedParts parts(rings, adjacency);
    const BoundaryLoops ringLoops(adjacency, parts);

    const bool locked = frozen_->lockedEdges != nullptr && !frozen_->lockedEdges->empty();
    std::vector<int> ring;
    std::vector<char> mark;
    std::vector<int> run;
    for (int r = 0; r < ringLoops.count(); ++r) {
        const auto loop = ringLoops.points().row(r);
        ring.assign(loop.rbegin(), loop.rend());
        const std::size_t n = ring.size();
        std::size_t last = n;
        int markCount = 0;
        if (locked) {
            mark.assign(n, 0);
            for (std::size_t j = 0; j < n; ++j) {
                if (marked(ring[j])) {
                    mark[j] = 1;
                    last = j;
                    ++markCount;
                }
            }
        }
        if (markCount == 0) {
            if (n >= 3) {
                addPolygon(ring);
            }
            continue;
        }
        if (markCount == 1) {
            // One mark: the whole ring from the mark.
            std::rotate(ring.begin(), ring.begin() + static_cast<std::ptrdiff_t>(last), ring.end());
            if (n >= 3) {
                addPolygon(ring);
            }
            continue;
        }
        // From the last mark round to it again, cut at every mark.
        run.clear();
        for (std::size_t j = last; j <= last + n; ++j) {
            const std::size_t k = j % n;
            run.push_back(ring[k]);
            if (mark[k]) {
                if (run.size() >= 3) {
                    addPolygon(run);
                }
                run.assign(1, ring[k]);
            }
        }
    }
}

// A ring point on a locked border edge: from a point with a locked edge or
// from a locked edge of the triangle mesh.
bool HoleClosing::marked(int point) const
{
    if (point < 0 || static_cast<std::size_t>(point) >= points_.source.size()) {
        return false;
    }
    const PointSource source = points_.source[static_cast<std::size_t>(point)];
    const BitArray& locked = *frozen_->lockedEdges;
    if (source.type == 1) {
        for (const int e : original_.adjacency().pointEdges(source.id)) {
            if (locked.testBounded(static_cast<std::size_t>(e))) {
                return true;
            }
        }
        return false;
    }
    return source.type == 2 && locked.testBounded(static_cast<std::size_t>(source.id));
}

void HoleClosing::addPolygon(std::span<const int> polygon)
{
    polygonCorners_.insert(polygonCorners_.end(), polygon.begin(), polygon.end());
    polygonEnds_.push_back(static_cast<int>(polygonCorners_.size()));
}

// d) Old faces, then the new polygons; joined corners take their target.
void HoleClosing::apply()
{
    if (polygonEnds_.empty() && joinedCorners_.empty()) {
        return;
    }
    const Topology& topology = output_.topology();
    std::vector<int> faceStart{0};
    for (int f = 0; f < topology.faceCount(); ++f) {
        faceStart.push_back(topology.faceEnd(f));
    }
    std::vector<int> corners = topology.cornerVertex();
    const int offset = static_cast<int>(corners.size());
    for (const int end : polygonEnds_) {
        faceStart.push_back(offset + end);
    }
    corners.insert(corners.end(), polygonCorners_.begin(), polygonCorners_.end());
    for (std::size_t k = 0; k < joinedCorners_.size(); ++k) {
        corners[static_cast<std::size_t>(joinedCorners_[k])] = joinTargets_[k];
    }
    Modification change;
    change.faceStart = std::move(faceStart);
    change.cornerVertex = std::move(corners);
    change.faceMode = FaceMode::Polygons;
    change.points.state = 1;
    MeshModifier::apply(output_, change);
}

} // namespace areko::remesher::extract
