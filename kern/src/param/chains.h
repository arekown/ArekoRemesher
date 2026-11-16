#pragma once

#include "field/cross_field.h"
#include "mesh/bit_array.h"
#include "mesh/csr.h"

#include <cstdint>
#include <vector>

namespace areko::remesher::param {

// Chains of curve edges on the cut copy of a part: along each of them one
// coordinate of the parameterisation is to stay constant.
struct CurveChains {
    // Per edge of the cut copy the kind of the curve edge it lies on
    // (field::curve_edge), none for all others and for the curve edges that
    // get no chain (CurveChainBuilder).
    std::vector<std::uint8_t> edgeKind;
    // Chain -> edges of the cut copy. A chain keeps one direction label from
    // end to end; a run of curve edges is split where the label changes.
    mesh::Csr edges;
    // Per chain 1 when it runs along the first combed direction (its second
    // coordinate stays constant), 0 along the second (its first stays).
    std::vector<std::int8_t> coordinate;
    // Points of the part mesh where a run of curve edges was split, in the
    // order found. The extraction starts from them as from singular points.
    std::vector<int> turnPoints;

    [[nodiscard]] int count() const noexcept { return edges.rows(); }
};

// Chains of the cut graph of one part: the cut edges as runs between nodes,
// each with the points of the cut copy on its two sides and the quarter turns
// between them. The parameterisation ties the two sides of every chain
// together.
struct Chains {
    // Singular points the parameterisation works with: those of the cross
    // field. Without any, one node stands in with index 0 (`substituted`).
    field::Singularities singularities;
    bool substituted = false;

    // Nodes of the cut graph: the ends of the cuts, then their branch points,
    // then the singular points and the border points in the middle of a run.
    // Chains begin and end only at nodes (a closed loop without any begins at
    // the first node, or runs free of them).
    std::vector<int> nodes;
    // Per node whether it is a singular point (bit per entry of `nodes`).
    mesh::BitArray singularNodes;

    // Chain -> edges of the part mesh, and chain -> its points: a chain that
    // ends at a node has one point more than edges, a closed chain repeats
    // its first point.
    mesh::Csr edges;
    mesh::Csr points;

    // Per chain the points of the cut copy on its two sides, one entry per
    // chain point whose chain edge has two corners (ChainSides). Side A is
    // the side of the triangle of the first corner of the first chain edge.
    // Both are empty when the sides could not be followed along the chain (a
    // border edge behind the first edge) and for a chain of one edge without
    // exactly two corners.
    std::vector<std::vector<int>> sideA;
    std::vector<std::vector<int>> sideB;
    // Per chain the quarter turns between the combed crosses of the two
    // sides, seen from side A (0 for a chain without sides).
    std::vector<float> rotation;

    // Chains of one edge that hang on a singular point around which the
    // crosses make a whole turn, with an end that the cut does not split
    // (bit per chain): they get no transition.
    mesh::BitArray fullTurnChains;
    // Chains of one curve check cut without a turn (bit per chain): they
    // only tie the coordinate across the curve.
    mesh::BitArray curveCheckChains;

    // Chains along the disc cuts, ascending; the coordinates on their two
    // sides differ by an unknown offset. Per chain the number of its first
    // offset variable (the second follows it), 0 for a chain without.
    std::vector<int> genusChains;
    std::vector<int> offsetVariable;

    CurveChains curves;

    [[nodiscard]] int count() const noexcept { return edges.rows(); }
};

} // namespace areko::remesher::param
