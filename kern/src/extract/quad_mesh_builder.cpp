#include "quad_mesh_builder.h"

#include "face_closing.h"

#include "mesh/boundary_loops.h"
#include "parts/part.h"

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

// Calls visit(node) for every node, triangle by triangle ascending, in the
// list order of each triangle (newest first).
template <class Visit>
void forEachNode(const TraceNodes& nodes, Visit visit)
{
    for (int t = 0; t < nodes.triangleCount(); ++t) {
        nodes.findIn(t, [&](int n) {
            visit(n);
            return false;
        });
    }
}

} // namespace

RawMesh QuadMeshBuilder::build(const parts::Part& part, const Trace& trace)
{
    const Mesh& triangles = part.mesh;
    const TraceNodes& nodes = trace.nodes;
    const int count = nodes.idCount();
    const auto size = static_cast<std::size_t>(count);
    RawMesh raw;
    OutputPoints& points = raw.points;

    // Points of the triangle mesh with a node on them.
    raw.originalToOutput.assign(static_cast<std::size_t>(triangles.pointCount()), -1);
    forEachNode(nodes, [&](int n) {
        const TraceNode& node = nodes[n];
        const int target = FaceClosing::number(node);
        if (node.kind == node_kind::point && target >= 0) {
            raw.originalToOutput[static_cast<std::size_t>(node.reference)] = target;
        }
    });

    // Kept border points without a node get output points of their own.
    std::vector<Vec3f> extraPoints;
    std::vector<PointSource> extraSources;
    if (const std::vector<FreezeMark>* marks = part.frozenBorderMarks()) {
        const BoundaryLoops loops(triangles);
        for (const int v : loops.points().entries) {
            int& out = raw.originalToOutput[static_cast<std::size_t>(v)];
            if ((*marks)[static_cast<std::size_t>(v)].value == 1 && out == -1) {
                out = count + static_cast<int>(extraPoints.size());
                extraPoints.push_back(triangles.points()[static_cast<std::size_t>(v)]);
                extraSources.push_back({node_kind::point, v});
            }
        }
    }

    // Merged nodes: their numbers in the faces go to their targets.
    points.weight.assign(size, 0.0f);
    std::vector<int> renumber(size, -1);
    for (int n = 0; n < nodes.count(); ++n) {
        const TraceNode& node = nodes[n];
        if (node.id >= 0 && node.mergeTarget >= 0) {
            renumber[static_cast<std::size_t>(node.id)] = node.mergeTarget;
            points.weight[static_cast<std::size_t>(node.mergeTarget)] = 0.3f;
        }
    }
    std::vector<int> faceStart = trace.polygons.starts;
    std::vector<int> corners = trace.polygons.entries;
    for (int& e : corners) {
        const int target = renumber[static_cast<std::size_t>(e)];
        if (target != -1) {
            e = target;
        }
    }

    std::vector<Vec3f> positions(size, Vec3f{});
    points.source.assign(size, PointSource{0, 0});
    points.triangle.assign(size, 0);
    for (std::vector<int>& list : points.neighbours) {
        list.assign(size, 0);
    }
    points.reachedBorder = BitArray(size);
    points.onBorderEdge = BitArray(size);
    const double near = static_cast<double>(diagonal(triangles.bounds())) * 5e-6;
    const double nearSq = near * near;
    forEachNode(nodes, [&](int n) {
        TraceNode node = nodes[n];
        if (node.id < 0) {
            return;
        }
        const auto id = static_cast<std::size_t>(node.id);
        snapToCorner(triangles, nearSq, node);
        // The output point takes the place of the node after it moved onto a corner.
        positions[id] = node.position;
        points.source[id] = node.kind == node_kind::free ? PointSource{4, node.triangle}
                                                          : PointSource{node.kind, node.reference};
        points.triangle[id] = node.triangle;
        for (std::size_t s = 0; s < 4; ++s) {
            const int other = node.links[s].node;
            points.neighbours[s][id] = other < 0 ? -1 : FaceClosing::number(nodes[other]);
        }
        if (node.seed()) {
            const bool onCurve = node.kind == node_kind::point &&
                                 trace.curves.points.testBounded(static_cast<std::size_t>(node.reference));
            points.weight[id] = onCurve ? 0.0f : 1.0f;
        }
        if ((node.flags & (node_flag::border | node_flag::ranIntoBorder)) != 0) {
            points.reachedBorder.set(id);
        }
        if ((node.flags & node_flag::border) != 0) {
            points.onBorderEdge.set(id);
        }
    });

    positions.insert(positions.end(), extraPoints.begin(), extraPoints.end());
    points.source.insert(points.source.end(), extraSources.begin(), extraSources.end());
    points.weight.resize(positions.size(), 0.0f);
    raw.mesh.swapPoints(positions);
    raw.mesh.editTopology().swapPolygons(faceStart, corners);
    raw.mesh.topologyChanged();
    return raw;
}

void QuadMeshBuilder::snapToCorner(const Mesh& triangles, double toleranceSq, TraceNode& node)
{
    if (node.triangle < 0) {
        return;
    }
    const Vec3f from = node.position;
    const Topology& topology = triangles.topology();
    double limit = toleranceSq;
    for (int c = topology.faceBegin(node.triangle); c < topology.faceEnd(node.triangle); ++c) {
        const int v = topology.vertex(c);
        const Vec3f corner = triangles.points()[static_cast<std::size_t>(v)];
        const Vec3f d = from - corner;
        const double d2 = static_cast<double>(dot(d, d));
        if (limit > d2) {
            node.kind = node_kind::point;
            node.reference = v;
            node.position = corner;
            limit = d2;
        }
    }
}

} // namespace areko::remesher::extract
