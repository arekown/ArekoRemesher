#include "curve_marks.h"

#include "areko_remesher.h"
#include "field/curve_edges.h"
#include "mesh/edge_chains.h"
#include "parts/part.h"

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

CurveMarks CurveMarks::mark(const parts::Part& part)
{
    const Adjacency& adjacency = part.mesh.adjacency();
    const auto edgeTotal = static_cast<std::size_t>(adjacency.edgeCount());
    CurveMarks out;
    out.points = BitArray(static_cast<std::size_t>(part.mesh.pointCount()));
    out.edges = BitArray(edgeTotal);
    for (const curves::GuideCurve& curve : part.curves) {
        if (!curve.active || curve.kind != ARM_CURVE_CREASE) {
            continue;
        }
        for (const int v : curve.cutVertices) {
            out.points.setBounded(static_cast<std::size_t>(v));
        }
        field::CurveEdges::forEachChainEdge(adjacency, curve, [&](int e) {
            if (!out.edges.test(static_cast<std::size_t>(e))) {
                out.edges.set(static_cast<std::size_t>(e));
                ++out.edgeCount;
            }
        });
    }
    out.chainPoints = EdgeChains::build(adjacency, out.edges, &out.chainEdges);

    out.edgeChain.assign(edgeTotal, -1);
    out.edgePlace.assign(edgeTotal, 0);
    out.alongFirstEnd = BitArray(edgeTotal);
    for (int chain = 0; chain < out.chainEdges.rows(); ++chain) {
        const auto edges = out.chainEdges.row(chain);
        const auto points = out.chainPoints.row(chain);
        for (std::size_t s = 0; s < edges.size(); ++s) {
            const auto e = static_cast<std::size_t>(edges[s]);
            out.edgeChain[e] = chain;
            out.edgePlace[e] = static_cast<int>(s);
            if (points[s] == adjacency.edgeVertices(edges[s])[0]) {
                out.alongFirstEnd.set(e);
            }
        }
    }
    return out;
}

} // namespace areko::remesher::extract
