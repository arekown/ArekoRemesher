#pragma once

#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Welds open seams: border points with exactly equal coordinates become one
// point (no tolerance). Candidates are the end points of the border edges,
// border group by border group, each point once at its first appearance;
// points with a NaN coordinate never weld. Candidates are sorted by x, y, z;
// neighbours in this order with equal coordinates join the first point of
// their run. Sharp edges keep their place on the welded mesh; with `asCrease`
// the welded former border edges become sharp edges as well.
class SeamWelder {
public:
    // Returns the number of points that were moved onto another one.
    static int apply(mesh::Mesh& mesh, bool asCrease);
};

} // namespace areko::remesher::prep
