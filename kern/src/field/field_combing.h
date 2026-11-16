#pragma once

#include "cross_field.h"
#include "cut_graph.h"

#include "mesh/mesh.h"

#include <vector>

namespace areko::remesher::field {

// Combs the cross field over the cut mesh: every triangle gets the number of
// quarter turns that lines its cross up with its neighbours.
//  - from the smallest triangle not yet visited the turns spread breadth
//    first across every edge that is not cut, taking the quarter turns of the
//    cross field along,
//  - the mirrored half takes 2 minus the turns of its mirror triangles,
//  - the combed directions are the first frame axis turned about the normal
//    by the field angle plus the quarter turns, and that direction turned by
//    another quarter, both written in the frame of the triangle.
class FieldCombing {
public:
    // Writes rotation, combedU and combedV of `out`. `faceMirror`: mirror
    // triangle per triangle when the part is symmetric, else null.
    static void comb(const mesh::Mesh& mesh, const CrossField& field, const std::vector<int>* faceMirror, CutGraph& out);

    // Writes residualJump of `out` from the rotation written by comb(): per
    // edge with exactly two corners the turns of its second triangle minus
    // those of its first plus the quarter turns of the cross field, else 0.
    static void residualJumps(const mesh::Adjacency& adjacency, const CrossField& field, CutGraph& out);
};

} // namespace areko::remesher::field
