#pragma once

#include "mesh.h"

namespace areko::remesher::mesh {

// Removes degenerate faces and corners of a mesh:
//  - faces with fewer than three corners and quads whose corners 0 and 2 or
//    1 and 3 name the same point vanish,
//  - a corner without point (-1) is left out, as is a corner whose point
//    equals the point of the corner before it in the old face (also when that
//    one was -1 or left out); the first corner is compared with the last
//    corner of the face that names a point,
//  - a face left with fewer than three corners vanishes; the others keep
//    their order.
// Only a real change touches the mesh; unused points stay. How the change
// reaches the maps:
//  - keepMaps: the kept old corners of every face go through FaceRebuild:
//    the points and their maps stay, faces, corners and half edges follow
//    their origin (the hard edges follow the half edges), and
//    with a face map the face mode stays, otherwise the faces are polygons.
//    The preparation cleans welded border loops this way, and it is the way
//    for the cleanup at the end of the extraction, which keeps hard edges
//    and point maps;
//  - otherwise without face origin the faces are set directly as polygons
//    and nothing is renumbered (point maps stay; face, edge and half edge
//    maps go), as right after the raw mesh is built;
//  - otherwise the face origin, and with the corner origin the old corners
//    that were not left out (ascending, brought to the new corner count) as
//    corner and half edge origin, go to the mesh change.
class DegenerateCleanup {
public:
    struct Switches {
        // Hand the face origin to the change when the mesh has a face map
        // (without keepMaps).
        bool faceOrigin = false;
        // Hand the corner origin as half edge origin when the mesh has an
        // edge or half edge map (without keepMaps).
        bool cornerOrigin = false;
        bool keepMaps = false;
    };

    struct Result {
        int removedFaces = 0;
    };

    static Result run(Mesh& mesh, const Switches& switches);
};

} // namespace areko::remesher::mesh
