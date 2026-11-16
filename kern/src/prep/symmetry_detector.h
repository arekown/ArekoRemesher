#pragma once

#include "symmetry_cache.h"

#include "mesh/mesh.h"

#include <optional>

namespace areko::remesher::prep {

// Finds the mirror symmetry of a triangle mesh from its topology:
//  - start edges lie almost parallel to the middle plane of the bounding box
//    along x, y or z (a locked plane of `previous` allows only its axis), in a
//    round with a strict and one without a direction limit,
//  - from a start edge the line where the plane cuts the mesh is followed in
//    both directions, along edges and through faces; edges once on a line or
//    crossed are not taken again; lines with too many changes between edges
//    and faces are dropped, and so is a line that leaves the plane and turns
//    by more than 10 degrees in the tangent plane of the split normals,
//  - from the faces of one side at the line a breadth first search pairs
//    faces with mirror faces and with them edges and points; the start edge
//    with the most paired faces wins, a full pairing ends the search,
//  - more than 70 percent paired faces make the mesh symmetric: side 1 are
//    the corners of the searched faces, the normal is the mean direction
//    between paired points (taken as the axis when a component exceeds
//    0.995), and sides are turned so that side 1 lies towards the normal; the
//    point is the mean of all paired points (a locked plane stays as it is),
//  - when not every face is paired, points without partner are paired across
//    the plane by distance; edge and face pairs stay those of the search.
// Meshes in other modes and meshes without symmetry give no result.
class SymmetryDetector {
public:
    [[nodiscard]] static std::optional<SymmetryCache> run(const mesh::Mesh& mesh, const SymmetryCache* previous);
};

} // namespace areko::remesher::prep
