#pragma once

#include "cut_graph.h"

#include "parts/part.h"
#include "prep/symmetry_cache.h"

namespace areko::remesher::field {

// Collects the edges of a part mesh that the quads are to follow, in this
// order (an edge may appear several times; its kind is the last one written):
//  - the edges on the mirror plane whose first end point belongs to a piece
//    of the mesh, when the part is evaluated as symmetric,
//  - all border edges, when borders are followed,
//  - the edges along the cut chain of every active curve that is cut into
//    the mesh, with the kind of the curve.
class CurveEdges {
public:
    struct Options {
        bool followBorders = false;
    };

    // `symmetric`: the mirror pairs of a symmetric part whose evaluation is
    // valid, else null. Writes curveEdges and curveEdgeKind of `out`.
    static void collect(const parts::Part& part, prep::SymmetryCache* symmetric, const Options& options, CutGraph& out);

    // Calls `visit(e)` for the edge e between every two consecutive points of
    // the cut chain of `curve` that an edge joins, in chain order. A point
    // twice in a row joins nothing (no edge of its own is taken for it).
    template <class Visit>
    static void forEachChainEdge(const mesh::Adjacency& adjacency, const curves::GuideCurve& curve, Visit&& visit)
    {
        const std::vector<int>& chain = curve.cutVertices;
        for (std::size_t i = 0; i + 1 < chain.size(); ++i) {
            const int e = adjacency.edgeBetween(chain[i], chain[i + 1]);
            if (e != -1) {
                visit(e);
            }
        }
    }
};

} // namespace areko::remesher::field
