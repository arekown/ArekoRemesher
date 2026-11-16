#pragma once

#include "attribute_map.h"
#include "bit_array.h"
#include "csr.h"
#include "index_map.h"

#include <optional>

namespace areko::remesher::mesh {

class Mesh;

// The sharp edge map exists either per edge or per half edge (one bit per
// corner). Edge -> corners sets every corner of a set edge; corners -> edge
// sets an edge when one of its corners is set.
class SharpEdges {
public:
    [[nodiscard]] static BoolMap* edgeForm(Mesh& mesh);
    // Edge form -> half edge form; true when a conversion took place.
    static bool toHalfEdgeForm(Mesh& mesh);
    // Half edge form -> edge form on the current topology; true when converted.
    static bool toEdgeForm(Mesh& mesh);
    // Sets the edge form, removes the half edge form.
    static void setEdgeBits(Mesh& mesh, BitArray bits);

    // Read without changing the mesh, nothing when there is no map:
    //  - one bit per corner: the half edge form as stored, else the edge form,
    //  - one bit per edge of the current topology: the edge form, else the
    //    half edge form.
    [[nodiscard]] static std::optional<BitArray> halfEdgeBits(const Mesh& mesh);
    [[nodiscard]] static std::optional<BitArray> edgeBits(const Mesh& mesh);
    // Start for adding sharp edges: the existing ones (a half edge form is
    // turned into the edge form first) or none, one bit per edge of the current
    // topology. The result goes back with setEdgeBits.
    [[nodiscard]] static BitArray extendableEdgeBits(Mesh& mesh);

    // Lays point chains of the old mesh onto the current one. Points keep their
    // numbers when the point state is 1 or 2 (only old points count), otherwise
    // they follow the inverted new -> old point list; without that list the map
    // stays as it is. Consecutive points a, b mark the edge (a, b) or, without
    // it, the two edges over the common neighbour c with the smallest
    // |a-c| + |c-b| (first one on ties).
    static void fromChains(Mesh& mesh, const Csr& chains, const IndexMap& points, int oldPointCount);
};

} // namespace areko::remesher::mesh
