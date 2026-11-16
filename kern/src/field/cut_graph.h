#pragma once

#include "cut_copy.h"

#include "areko_remesher.h"
#include "mesh/adjacency.h"
#include "mesh/bit_array.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace areko::remesher::field {

// Kinds of curve edges: edges the quads are to follow. An edge along a curve
// that is cut into the mesh carries the kind of that curve.
namespace curve_edge {
inline constexpr std::uint8_t none = 0;
inline constexpr std::uint8_t edgeLoop = ARM_CURVE_EDGE_LOOP;
inline constexpr std::uint8_t crease = ARM_CURVE_CREASE; // also border edges and edges on the mirror plane
static_assert(edgeLoop != none && crease != none);
} // namespace curve_edge

// Calls visit(e) for the edges of the triangles set in `triangles`, in the
// edge numbers of `adjacency`: triangles ascending, their three corners in
// order, corners without an edge left out. An edge between two such
// triangles comes twice.
template <class Visit>
void forEachEdgeOfTriangles(const mesh::Adjacency& adjacency, const mesh::BitArray& triangles, Visit visit)
{
    triangles.forEachSet([&](std::size_t t) {
        for (int i = 0; i < 3; ++i) {
            const int e = adjacency.cornerEdge(3 * static_cast<int>(t) + i);
            if (e >= 0) {
                visit(e);
            }
        }
    });
}

// Cut graph of one part: the edges along which the part is cut open into a
// disc, and the cross field combed over that disc. Later steps only read it.
struct CutGraph {
    // All cut edges of the part mesh (bit per edge, ascending list) and the
    // mesh cut open along them.
    mesh::BitArray cutEdgeBits;
    std::vector<int> cutEdges;
    CutCopy copy;
    // Disc cuts: the cuts among them that were added last, where the copy cut
    // open along tree, plane chains and stars was not yet a disc (found by
    // DiskCut). Every pass finds its own; the bit array is empty without one.
    std::vector<int> diskCuts;
    mesh::BitArray diskCutBits;
    // Curve check cuts: edges the curve check (CurveAlignment::check) adds to
    // the cuts where two curves at a point would get the same direction, and
    // the points it adds as sources of the tree. They are kept over the
    // passes; the bit array is empty until the first one.
    std::vector<int> curveCheckCuts;
    mesh::BitArray curveCheckCutBits;
    std::vector<int> curveCheckSources;

    // Per triangle the number of quarter turns that brings its cross in line
    // with its neighbours across every edge that is not cut (0 .. 3; in the
    // mirrored half 2 minus the value of the mirror triangle, so -1 .. 2), and
    // the two combed directions in the frame of the triangle.
    std::vector<int> rotation;
    std::vector<std::array<float, 2>> combedU;
    std::vector<std::array<float, 2>> combedV;
    // Per edge the quarter turns left between the combed crosses of its two
    // triangles (0 for edges without exactly two corners).
    std::vector<int> residualJump;

    // Curve edges (an edge may be listed more than once), their kind per edge
    // (curve_edge), and per corner whether the curve edge starting there runs
    // along the first combed direction rather than the second.
    std::vector<int> curveEdges;
    std::vector<std::uint8_t> curveEdgeKind;
    mesh::BitArray curveFollowsU;

    // Triangles of the mirrored half and how many they are.
    mesh::BitArray mirrorHalf;
    int mirrorHalfCount = 0;
};

} // namespace areko::remesher::field
