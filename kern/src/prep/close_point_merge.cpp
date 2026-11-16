#include "close_point_merge.h"

#include "mesh/edge_groups.h"
#include "mesh/face_rebuild.h"

#include <utility>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

int ClosePointMerge::apply(Mesh& mesh, float nearSquared, float smallSquared)
{
    const Adjacency& adjacency = mesh.adjacency();
    const auto& points = mesh.points();
    const bool triangles = mesh.topology().mode() == FaceMode::Triangles && smallSquared > 0.0f;
    float limit = 0.0f;
    const std::vector<float>* doubleAreas = nullptr;
    if (triangles) {
        const TriangleGeometry& geometry = mesh.triangleGeometry();
        limit = geometry.smallAreaLimit();
        doubleAreas = &geometry.doubleAreas();
    }

    std::vector<int> candidates;
    for (int e = 0; e < adjacency.edgeCount(); ++e) {
        const auto ends = adjacency.edgeVertices(e);
        const Vec3f a = points[static_cast<std::size_t>(ends[0])];
        const Vec3f b = points[static_cast<std::size_t>(ends[1])];
        if (!triangles && nearSquared == 0.0f) {
            if (a == b) {
                candidates.push_back(e);
            }
            continue;
        }
        const Vec3f d = a - b;
        const float l2 = dot(d, d);
        if (l2 <= nearSquared) {
            candidates.push_back(e);
        } else if (triangles && l2 <= smallSquared) {
            // Every triangle at the edge must be small (also at edges with more than two corners).
            bool small = true;
            for (const int c : adjacency.edgeCorners(e)) {
                small = small && !((*doubleAreas)[static_cast<std::size_t>(adjacency.cornerFace(c))] > limit);
            }
            if (small) {
                candidates.push_back(e);
            }
        }
    }
    if (candidates.empty()) {
        return 0;
    }

    const Csr groups = EdgeGroups::byConnectivity(adjacency, candidates);
    std::vector<int> target(static_cast<std::size_t>(mesh.pointCount()), -1);
    std::vector<int> sources;
    std::vector<int> unique;
    for (int g = 0; g < groups.rows(); ++g) {
        unique.clear();
        for (const int e : groups.row(g)) {
            for (const int x : adjacency.edgeVertices(e)) {
                if (target[static_cast<std::size_t>(x)] != -2) {
                    target[static_cast<std::size_t>(x)] = -2;
                    unique.push_back(x);
                }
            }
        }
        if (unique.size() < 2) {
            continue;
        }
        for (std::size_t i = 0; i < unique.size(); ++i) {
            if (adjacency.firstBoundaryEdgeAt(unique[i]) >= 0) {
                std::swap(unique[0], unique[i]);
                break;
            }
        }
        for (std::size_t i = 1; i < unique.size(); ++i) {
            sources.push_back(unique[i]);
            target[static_cast<std::size_t>(unique[i])] = unique[0];
        }
    }
    if (sources.empty()) {
        return 0;
    }

    const Topology& topology = mesh.topology();
    std::vector<int> corners = topology.cornerVertex();
    BitArray touched(static_cast<std::size_t>(topology.faceCount()));
    for (const int q : sources) {
        for (const int c : adjacency.pointCorners(q)) {
            corners[static_cast<std::size_t>(c)] = target[static_cast<std::size_t>(q)];
            touched.set(static_cast<std::size_t>(adjacency.cornerFace(c)));
        }
    }
    FaceRebuild rebuild(mesh, std::move(corners));
    std::vector<int> kept;
    for (int f = 0; f < topology.faceCount(); ++f) {
        if (!touched.test(static_cast<std::size_t>(f))) {
            rebuild.keepAll(f);
            continue;
        }
        kept.clear();
        const int b = topology.faceBegin(f);
        const int n = topology.faceSize(f);
        for (int i = 0; i < n; ++i) {
            if (rebuild.point(b + i) != rebuild.point(b + (i + 1) % n)) {
                kept.push_back(b + i);
            }
        }
        rebuild.keep(f, kept);
    }
    rebuild.apply(mesh, topology.mode());
    return static_cast<int>(sources.size());
}

} // namespace areko::remesher::prep
