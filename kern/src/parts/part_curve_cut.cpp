#include "part_curve_cut.h"

#include "curves/curve_projector.h"
#include "curves/segment_cutter.h"
#include "curves/segment_preparation.h"
#include "prep/edge_splitter.h"
#include "prep/sliver_fixer.h"

#include <vector>

namespace areko::remesher::parts {

using namespace areko::remesher::mesh;
using curves::CurveSegment;
using curves::GuideCurve;

void PartCurveCut::run(Part& part, bool mirrorAxis)
{
    cut(part);
    densify(part, mirrorAxis);
}

void PartCurveCut::cut(Part& part)
{
    std::vector<const CurveSegment*> pieces;
    for (GuideCurve& curve : part.curves) {
        if (!curves::cutsIntoMesh(curve.kind)) {
            continue;
        }
        if (!curve.projected) {
            curve.segments = curves::CurveProjector::project(part.mesh, curve.source);
            curve.projected = true;
        }
        for (const CurveSegment& piece : curve.segments) {
            pieces.push_back(&piece);
        }
    }
    if (pieces.empty()) {
        return;
    }
    const std::vector<CurveSegment> prepared =
        curves::SegmentPreparation::prepare(part.mesh, pieces, curves::SegmentPreparation::Options{});
    const curves::SegmentCutter::Result result =
        curves::SegmentCutter::cut(part.mesh, prepared, curves::SegmentCutter::Options{});
    if (result.changed) {
        part.mirror.invalidate();
    }
    std::size_t s = 0;
    for (GuideCurve& curve : part.curves) {
        if (!curve.active || !curves::cutsIntoMesh(curve.kind)) {
            continue;
        }
        curve.cutVertices.clear();
        curve.cutPoints.clear();
        for (std::size_t k = 0; k < curve.segments.size() && s < result.chains.size(); ++k, ++s) {
            for (const int v : result.chains[s]) {
                if (v < 0) {
                    continue;
                }
                curve.cutVertices.push_back(v);
                curve.cutPoints.push_back(part.mesh.points()[static_cast<std::size_t>(v)]);
            }
        }
    }
}

void PartCurveCut::densify(Part& part, bool mirrorAxis)
{
    Mesh& mesh = part.mesh;
    const Adjacency& adjacency = mesh.adjacency();
    const int pointCount = mesh.pointCount();
    BitArray featurePoints(static_cast<std::size_t>(pointCount));
    BitArray featureEdges(static_cast<std::size_t>(adjacency.edgeCount()));
    for (const GuideCurve& curve : part.curves) {
        if (!curves::cutsIntoMesh(curve.kind)) {
            continue;
        }
        const std::vector<int>& chain = curve.cutVertices;
        for (std::size_t i = 0; i + 1 < chain.size(); ++i) {
            const int a = chain[i];
            const int b = chain[i + 1];
            if (a == -1 || b == -1) {
                continue;
            }
            featurePoints.setBounded(static_cast<std::size_t>(a));
            featurePoints.setBounded(static_cast<std::size_t>(b));
            const int e = adjacency.edgeBetween(a, b);
            if (e != -1) {
                featureEdges.set(static_cast<std::size_t>(e));
            }
        }
    }
    for (const int e : adjacency.boundaryEdges()) {
        featureEdges.set(static_cast<std::size_t>(e));
    }
    for (const int v : adjacency.boundaryPoints()) {
        featurePoints.set(static_cast<std::size_t>(v));
    }
    prep::SymmetryCache* mirror = nullptr;
    if (mirrorAxis) {
        mirror = part.mirror.symmetric(mesh, true);
        if (mirror != nullptr) {
            for (const int e : mirror->planeEdges(adjacency)) {
                const auto ends = adjacency.edgeVertices(e);
                featurePoints.set(static_cast<std::size_t>(ends[0]));
                featurePoints.set(static_cast<std::size_t>(ends[1]));
                featureEdges.set(static_cast<std::size_t>(e));
            }
        }
    }
    std::vector<int> edges;
    for (int e = 0; e < adjacency.edgeCount(); ++e) {
        const auto ends = adjacency.edgeVertices(e);
        if (featurePoints.testBounded(static_cast<std::size_t>(ends[0])) &&
            featurePoints.testBounded(static_cast<std::size_t>(ends[1])) && !featureEdges.test(static_cast<std::size_t>(e))) {
            edges.push_back(e);
        }
    }
    BitArray faces;
    for (int f = 0; f < mesh.faceCount(); ++f) {
        bool enclosed = true;
        for (int c = 3 * f; c < 3 * f + 3 && enclosed; ++c) {
            const int e = adjacency.cornerEdge(c);
            enclosed = e >= 0 && featureEdges.test(static_cast<std::size_t>(e));
        }
        if (enclosed) {
            if (faces.empty()) {
                faces.resize(static_cast<std::size_t>(mesh.faceCount()));
            }
            faces.set(static_cast<std::size_t>(f));
        }
    }
    if (edges.empty() && faces.empty()) {
        return;
    }
    prep::EdgeSplitter::splitMarked(mesh, edges, faces.empty() ? nullptr : &faces, mirror, !mirrorAxis);
    part.mirror.invalidate();
    // Old points stay, new ones may move away from slivers.
    BitArray locked(static_cast<std::size_t>(mesh.pointCount()));
    for (int v = 0; v < pointCount; ++v) {
        locked.set(static_cast<std::size_t>(v));
    }
    prep::SliverFixer::run(mesh, locked);
}

} // namespace areko::remesher::parts
