#include "adjacency.h"

#include "bit_array.h"

namespace areko::remesher::mesh {

Adjacency::Adjacency(const Topology& topology, int pointCount) : topology_(&topology), pointCount_(pointCount)
{
    pointCorners_ = Csr::fromValues(topology.cornerVertex(), pointCount_);
    buildCornerFaces();
    buildCornerEdges();
    buildEdges();
    buildPointEdges();
}

void Adjacency::buildCornerFaces()
{
    const Topology& topology = *topology_;
    const int corners = topology.cornerCount();
    cornerFace_.assign(static_cast<std::size_t>(corners), -1);
    const int faces = topology.faceCount();
    for (int f = 0; f < faces; ++f) {
        const int end = std::min(topology.faceEnd(f), corners);
        for (int c = topology.faceBegin(f); c < end; ++c) {
            cornerFace_[static_cast<std::size_t>(c)] = f;
        }
    }
}

int Adjacency::next(int c) const noexcept
{
    switch (topology_->mode()) {
    case FaceMode::Triangles: {
        static constexpr int step[3] = {1, 1, -2};
        return c + step[c % 3];
    }
    case FaceMode::Quads:
        return (c - c % 4) + (c + 1) % 4;
    case FaceMode::Polygons: {
        const int f = cornerFace_[static_cast<std::size_t>(c)];
        const int n = c + 1;
        return n < topology_->faceEnd(f) ? n : topology_->faceBegin(f);
    }
    }
    return c;
}

int Adjacency::prev(int c) const noexcept
{
    switch (topology_->mode()) {
    case FaceMode::Triangles: {
        static constexpr int step[3] = {2, -1, -1};
        return c + step[c % 3];
    }
    case FaceMode::Quads:
        return (c - c % 4) + (c + 3) % 4;
    case FaceMode::Polygons: {
        const int f = cornerFace_[static_cast<std::size_t>(c)];
        return c > topology_->faceBegin(f) ? c - 1 : topology_->faceEnd(f) - 1;
    }
    }
    return c;
}

void Adjacency::buildCornerEdges()
{
    cornerEdge_.assign(cornerFace_.size(), -1);
    std::vector<int> mark(static_cast<std::size_t>(pointCount_), -1);
    int e = 0;
    for (int a = 0; a < pointCount_; ++a) {
        const int start = e;
        for (const int c : pointCorners_.row(a)) {
            if (cornerFace_[static_cast<std::size_t>(c)] < 0) {
                continue;
            }
            int b = vertex(next(c));
            if (b > a) {
                if (mark[static_cast<std::size_t>(b)] < start) {
                    mark[static_cast<std::size_t>(b)] = e++;
                }
                cornerEdge_[static_cast<std::size_t>(c)] = mark[static_cast<std::size_t>(b)];
            }
            const int p = prev(c);
            b = vertex(p);
            if (b > a) {
                if (mark[static_cast<std::size_t>(b)] < start) {
                    mark[static_cast<std::size_t>(b)] = e++;
                }
                cornerEdge_[static_cast<std::size_t>(p)] = mark[static_cast<std::size_t>(b)];
            }
        }
    }
    edgeCount_ = e;
}

void Adjacency::buildEdges()
{
    edgeVertices_.assign(2 * static_cast<std::size_t>(edgeCount_), -1);
    for (std::size_t c = 0; c < cornerEdge_.size(); ++c) {
        const int e = cornerEdge_[c];
        if (e < 0 || edgeVertices_[2 * static_cast<std::size_t>(e)] != -1) {
            continue;
        }
        const int a = vertex(static_cast<int>(c));
        const int b = vertex(next(static_cast<int>(c)));
        edgeVertices_[2 * static_cast<std::size_t>(e)] = a < b ? a : b;
        edgeVertices_[2 * static_cast<std::size_t>(e) + 1] = a < b ? b : a;
    }
    edgeCorners_ = Csr::fromValues(cornerEdge_, edgeCount_);
}

void Adjacency::buildPointEdges()
{
    // Counting sort over the end point list (2e, 2e+1); entry = e, ascending.
    pointEdges_.starts.assign(static_cast<std::size_t>(pointCount_) + 1, 0);
    for (const int v : edgeVertices_) {
        ++pointEdges_.starts[static_cast<std::size_t>(v) + 1];
    }
    for (int v = 0; v < pointCount_; ++v) {
        pointEdges_.starts[static_cast<std::size_t>(v) + 1] += pointEdges_.starts[static_cast<std::size_t>(v)];
    }
    pointEdges_.entries.resize(edgeVertices_.size());
    pointNeighbours_.resize(edgeVertices_.size());
    std::vector<int> fill(pointEdges_.starts.begin(), pointEdges_.starts.end() - 1);
    for (std::size_t i = 0; i < edgeVertices_.size(); ++i) {
        const int v = edgeVertices_[i];
        const int slot = fill[static_cast<std::size_t>(v)]++;
        pointEdges_.entries[static_cast<std::size_t>(slot)] = static_cast<int>(i >> 1);
        pointNeighbours_[static_cast<std::size_t>(slot)] = edgeVertices_[i ^ 1u];
    }
}

std::span<const int> Adjacency::pointNeighbours(int v) const
{
    const int begin = pointEdges_.starts[static_cast<std::size_t>(v)];
    return {pointNeighbours_.data() + begin, static_cast<std::size_t>(pointEdges_.rowSize(v))};
}

int Adjacency::otherEnd(int e, int v) const noexcept
{
    const auto ends = edgeVertices(e);
    return ends[0] == v ? ends[1] : ends[0];
}

int Adjacency::edgeBetween(int a, int b) const
{
    const auto edges = pointEdges(a);
    const auto others = pointNeighbours(a);
    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (others[i] == b) {
            return edges[i];
        }
    }
    return -1;
}

std::optional<std::array<int, 2>> Adjacency::innerEdgeFaces(int e) const
{
    const auto corners = edgeCorners(e);
    if (corners.size() != 2) {
        return std::nullopt;
    }
    return std::array<int, 2>{cornerFace(corners[0]), cornerFace(corners[1])};
}

int Adjacency::otherFace(int e, int f) const
{
    const auto faces = innerEdgeFaces(e);
    if (!faces) {
        return -1;
    }
    const auto [f0, f1] = *faces;
    return f0 == f ? f1 : (f1 == f ? f0 : -1);
}

int Adjacency::otherCorner(int e, int c) const
{
    const auto corners = edgeCorners(e);
    if (corners.size() != 2) {
        return -1;
    }
    return corners[0] == c ? corners[1] : (corners[1] == c ? corners[0] : -1);
}

int Adjacency::commonEdge(int f, int g) const
{
    if (f < 0 || f >= topology_->faceCount()) {
        return -1;
    }
    // The corners of f are the ones buildCornerFaces gave to it.
    const int end = std::min(topology_->faceEnd(f), static_cast<int>(cornerFace_.size()));
    for (int c = topology_->faceBegin(f); c < end; ++c) {
        const int e = cornerEdge(c);
        if (e < 0) {
            continue;
        }
        for (const int other : edgeCorners(e)) {
            if (cornerFace(other) == g) {
                return e;
            }
        }
    }
    return -1;
}

void Adjacency::buildBoundaryLists() const
{
    std::vector<int> boundary;
    nonManifoldEdges_.clear();
    for (int e = 0; e < edgeCount_; ++e) {
        const int n = edgeCorners_.rowSize(e);
        if (n < 2) {
            boundary.push_back(e);
        } else if (n > 2) {
            nonManifoldEdges_.push_back(e);
        }
    }
    boundaryEdges_ = std::move(boundary);
}

const std::vector<int>& Adjacency::boundaryEdges() const
{
    if (!boundaryEdges_) {
        buildBoundaryLists();
    }
    return *boundaryEdges_;
}

const std::vector<int>& Adjacency::nonManifoldEdges() const
{
    if (!boundaryEdges_) {
        buildBoundaryLists();
    }
    return nonManifoldEdges_;
}

const std::vector<int>& Adjacency::boundaryPoints() const
{
    if (!boundaryPoints_) {
        BitArray seen(static_cast<std::size_t>(pointCount_));
        std::vector<int> points;
        for (const int e : boundaryEdges()) {
            for (const int v : edgeVertices(e)) {
                if (!seen.test(static_cast<std::size_t>(v))) {
                    seen.set(static_cast<std::size_t>(v));
                    points.push_back(v);
                }
            }
        }
        boundaryPoints_ = std::move(points);
    }
    return *boundaryPoints_;
}

const std::vector<int>& Adjacency::isolatedPoints() const
{
    if (!isolatedPoints_) {
        std::vector<int> points;
        for (int v = 0; v < pointCount_; ++v) {
            if (pointCorners_.rowSize(v) == 0) {
                points.push_back(v);
            }
        }
        isolatedPoints_ = std::move(points);
    }
    return *isolatedPoints_;
}

int Adjacency::firstBoundaryEdgeAt(int v, int except) const
{
    for (const int e : pointEdges(v)) {
        if (isBoundaryEdge(e) && e != except) {
            return e;
        }
    }
    return -1;
}

int Adjacency::boundaryCornerAt(int v, bool outgoing) const
{
    for (const int c : pointCorners(v)) {
        if (cornerFace_[static_cast<std::size_t>(c)] < 0) {
            continue;
        }
        const int e = cornerEdge_[static_cast<std::size_t>(outgoing ? c : prev(c))];
        if (e >= 0 && isBoundaryEdge(e)) {
            return c;
        }
    }
    return -1;
}

const Adjacency::Counts& Adjacency::counts() const
{
    if (!counts_) {
        Counts n;
        const int faces = topology_->faceCount();
        for (int f = 0; f < faces; ++f) {
            const int size = topology_->faceSize(f);
            if (size == 4) {
                ++n.quads;
            } else if (size == 3) {
                ++n.triangles;
            } else {
                ++n.polygons;
            }
        }
        BitArray onBoundary(static_cast<std::size_t>(pointCount_));
        for (const int v : boundaryPoints()) {
            onBoundary.set(static_cast<std::size_t>(v));
        }
        for (int v = 0; v < pointCount_; ++v) {
            if (!onBoundary.test(static_cast<std::size_t>(v)) && pointCorners_.rowSize(v) > 0 &&
                pointEdges_.rowSize(v) != 4) {
                ++n.irregularInnerPoints;
            }
        }
        counts_ = n;
    }
    return *counts_;
}

} // namespace areko::remesher::mesh
