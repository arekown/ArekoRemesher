#pragma once

#include "raw_mesh.h"

#include "mesh/mesh.h"

namespace areko::remesher::extract {

// Face groups of the quad mesh from the face groups of the triangle mesh
// (nothing when the triangle mesh has none):
//  - every corner of a face takes a group from where its point comes from:
//    the triangle it lies in, or the point or edge of the triangle mesh when
//    all triangles there share one group; otherwise the triangle at that
//    point or edge that lies best in the direction of the face at the corner
//    (sum of the two sides of the corner),
//  - a face whose corners agree takes their group; the others stay open and
//    take, in up to 50 rounds, the group most of their neighbours across
//    edges that are no curve edges have (first one on ties); open faces
//    left over keep 0.
class GroupTransfer {
public:
    static void run(const mesh::Mesh& triangles, RawMesh& raw);

    // Without crease curves: every face takes the group of the first triangle.
    static void uniform(const mesh::Mesh& triangles, RawMesh& raw);
};

} // namespace areko::remesher::extract
