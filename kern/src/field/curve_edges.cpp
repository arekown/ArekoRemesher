#include "curve_edges.h"

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

void CurveEdges::collect(const parts::Part& part, prep::SymmetryCache* symmetric, const Options& options, CutGraph& out)
{
    const Mesh& mesh = part.mesh;
    const Adjacency& adjacency = mesh.adjacency();
    out.curveEdges.clear();
    out.curveEdgeKind.assign(static_cast<std::size_t>(adjacency.edgeCount()), curve_edge::none);
    auto add = [&](int e, std::uint8_t kind) {
        out.curveEdges.push_back(e);
        out.curveEdgeKind[static_cast<std::size_t>(e)] = kind;
    };
    if (symmetric != nullptr) {
        const ConnectedParts& pieces = mesh.connectedParts();
        for (const int e : symmetric->planeEdges(adjacency)) {
            if (pieces.partOfPoint(adjacency.edgeVertices(e)[0]) >= 0) {
                add(e, curve_edge::crease);
            }
        }
    }
    if (options.followBorders) {
        for (const int e : adjacency.boundaryEdges()) {
            add(e, curve_edge::crease);
        }
    }
    for (const curves::GuideCurve& curve : part.curves) {
        if (!curves::cutsIntoMesh(curve.kind) || !curve.active) {
            continue;
        }
        forEachChainEdge(adjacency, curve, [&](int e) { add(e, static_cast<std::uint8_t>(curve.kind)); });
    }
}

} // namespace areko::remesher::field
