#include "symmetry_cache.h"

#include <cstddef>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

int at(const std::vector<int>& values, int i)
{
    return i >= 0 && static_cast<std::size_t>(i) < values.size() ? values[static_cast<std::size_t>(i)] : -1;
}

} // namespace

SymmetryCache SymmetryCache::fromPoints(std::vector<int> points, BitArray pointSide,
                                        const std::optional<LockedPlane>& plane)
{
    SymmetryCache cache(Maps{std::move(points), std::move(pointSide), {}, {}}, LockedPlane{});
    cache.plane_ = plane;
    return cache;
}

void SymmetryCache::buildEdgeMirror(const Adjacency& adjacency)
{
    edgeMirrorBuilt_ = true;
    if (maps_.points.empty()) {
        edgeMirror_.clear();
        return;
    }
    const auto edges = static_cast<std::size_t>(adjacency.edgeCount());
    if (edgeMirror_.size() != edges) {
        edgeMirror_.assign(edges, -1);
    }
    for (int e = 0; e < adjacency.edgeCount(); ++e) {
        if (edgeMirror_[static_cast<std::size_t>(e)] >= 0) {
            continue;
        }
        const auto ends = adjacency.edgeVertices(e);
        const int a = at(maps_.points, ends[0]);
        const int b = at(maps_.points, ends[1]);
        if (a < 0 || b < 0 || a >= adjacency.pointCount()) {
            continue;
        }
        const int mirrored = adjacency.edgeBetween(a, b);
        if (mirrored < 0) {
            continue;
        }
        edgeMirror_[static_cast<std::size_t>(e)] = mirrored;
        edgeMirror_[static_cast<std::size_t>(mirrored)] = e;
    }
}

const std::vector<int>& SymmetryCache::edgeMirror(const Adjacency& adjacency)
{
    if (!edgeMirrorBuilt_) {
        buildEdgeMirror(adjacency);
    }
    return edgeMirror_;
}

const std::vector<int>& SymmetryCache::planeEdges(const Adjacency& adjacency)
{
    if (!planeEdges_) {
        planeEdges_.emplace();
        const std::vector<int>& mirror = edgeMirror(adjacency);
        for (std::size_t e = 0; e < mirror.size(); ++e) {
            if (mirror[e] == static_cast<int>(e)) {
                planeEdges_->push_back(static_cast<int>(e));
            }
        }
    }
    return *planeEdges_;
}

const std::vector<int>& SymmetryCache::faceMirror(const Mesh& mesh)
{
    if (mirroredFaces_ != -1) {
        return maps_.faces;
    }
    const auto faces = static_cast<std::size_t>(mesh.faceCount());
    if (maps_.points.empty()) {
        maps_.faces.clear();
        return maps_.faces;
    }
    maps_.faces.assign(faces, -1);
    if (mesh.topology().mode() != FaceMode::Triangles) {
        return maps_.faces;
    }
    const Adjacency& adjacency = mesh.adjacency();
    const Topology& topology = mesh.topology();
    mirroredFaces_ = 0;
    for (int f = 0; f < mesh.faceCount(); ++f) {
        if (maps_.faces[static_cast<std::size_t>(f)] >= 0) {
            continue;
        }
        const int a = at(maps_.points, topology.vertex(3 * f));
        const int b = at(maps_.points, topology.vertex(3 * f + 1));
        const int c = at(maps_.points, topology.vertex(3 * f + 2));
        if (a < 0 || b < 0 || c < 0 || a >= adjacency.pointCount()) {
            continue;
        }
        const int e = adjacency.edgeBetween(a, b);
        if (e < 0) {
            continue;
        }
        int g = -1;
        for (const int k : adjacency.edgeCorners(e)) {
            const int t = adjacency.cornerFace(k);
            if (topology.vertex(3 * t) == c || topology.vertex(3 * t + 1) == c || topology.vertex(3 * t + 2) == c) {
                g = t;
                break;
            }
        }
        if (g < 0) {
            continue;
        }
        // A target face that already has its pair is not taken again.
        if (g != f && maps_.faces[static_cast<std::size_t>(g)] >= 0) {
            continue;
        }
        maps_.faces[static_cast<std::size_t>(f)] = g;
        if (g == f) {
            mirroredFaces_ += 1;
        } else {
            maps_.faces[static_cast<std::size_t>(g)] = f;
            mirroredFaces_ += 2;
        }
    }
    return maps_.faces;
}

bool SymmetryCache::allFacesMirrored(const Mesh& mesh)
{
    (void)faceMirror(mesh);
    return mirroredFaces_ == mesh.faceCount();
}

int SymmetryCache::sideOfFace(const Topology& topology, int f) const
{
    if (maps_.points.empty()) {
        return -1;
    }
    int zero = 0;
    int one = 0;
    for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
        const int v = topology.vertex(c);
        if (at(maps_.points, v) == v) {
            continue;
        }
        if (maps_.pointSide.testBounded(static_cast<std::size_t>(v))) {
            ++one;
        } else {
            ++zero;
        }
    }
    if (zero == 0 && one > 0) {
        return 1;
    }
    if (zero > 0 && one == 0) {
        return 0;
    }
    if (zero == one && zero > 0) {
        return 2;
    }
    return -1;
}

} // namespace areko::remesher::prep
