#pragma once

#include "mesh/mesh.h"

namespace areko::remesher::merge {

// Splits the points of the result (a polygon mesh) along its sharp edges, so
// that the faces on either side of a sharp edge no longer share its points.
// Nothing happens without a set sharp edge; with one, the change is applied
// even when no point splits, so edge maps other than the sharp edges and UV
// edge maps are dropped then as well.
//  - cut edges are the sharp edges and the border edges; a point becomes a
//    candidate when its second cut edge is counted (sharp edges ascending,
//    then border edges ascending; the count per point is one byte),
//  - per candidate its faces join over every edge at it that is not sharp
//    (all faces of such an edge with the face of its first corner),
//  - the first fan (the one with the first corner of the point) keeps the
//    point, every further fan gets a new point at the same place, numbered on
//    from the point count in candidate and fan order; point maps copy the
//    point, faces and corners stay, a sharp edge stays sharp on both sides.
class CreaseSplit {
public:
    // Returns the number of new points.
    static int run(mesh::Mesh& mesh);
};

} // namespace areko::remesher::merge
