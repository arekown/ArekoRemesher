#pragma once

#include "chains.h"

#include "field/cross_field.h"
#include "field/cut_graph.h"
#include "mesh/mesh.h"

namespace areko::remesher::param {

// Builds the chains of curve edges of a part on its cut copy.
//  - around every point of a curve the direction label of one curve edge is
//    carried through the triangles to the next one, as in the curve check
//    (a cut edge turns it); a curve edge that arrives with the label of the
//    one before it is excluded, unless a curve check cut lies between them or
//    the sector between them is wider than a straight angle and bends back,
//  - the edges of the cut copy on the remaining curve edges, in the order of
//    the curve edge list, are sorted by the combed direction they follow and
//    chained, first those along the second direction, then those along the
//    first,
//  - every such run is split where the direction label of its edges changes
//    (seen in the triangle of the first corner of each edge, along the run);
//    the point of the split is noted. The last piece of a closed run joins
//    the first piece in front when both carry the same label,
//  - with a mirrored half, runs that begin in it are left out.
class CurveChainBuilder {
public:
    // `symmetry`: the part may be computed by halves.
    [[nodiscard]] static CurveChains build(const mesh::Mesh& mesh, const field::CrossField& field,
                                           const field::CutGraph& graph, bool symmetry);
};

} // namespace areko::remesher::param
