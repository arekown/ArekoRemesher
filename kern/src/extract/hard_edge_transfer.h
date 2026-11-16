#pragma once

#include "raw_mesh.h"

#include "curves/guide_curve.h"
#include "mesh/mesh.h"

#include <span>

namespace areko::remesher::extract {

// Sharp edges of the quad mesh from the sharp lines of the triangle mesh.
// Sharp lines are the active crease curves whose line flag has bit 4 or the
// low two bits equal to 1: the edges of the triangle mesh between their
// consecutive cut chain points. Without such an edge nothing happens. Every
// curve edge q of the quad mesh (ascending) is sharp when
//  - no end point comes from inside a triangle, and
//  - the first end point lies on an edge of the triangle mesh: that edge is
//    sharp; else the second one does: that edge is sharp; else both lie on
//    points a and b of the triangle mesh: among the edges at a that point the
//    same way as b - a (positive dot product) and within 45 degrees of it, the
//    one with the smallest angle is sharp – angles of sharp edges count 0.6
//    times, and a sharp edge within 6 degrees wins at once.
// The result is the sharp edge map of the quad mesh.
class HardEdgeTransfer {
public:
    static void run(const mesh::Mesh& triangles, std::span<const curves::GuideCurve> curves, RawMesh& raw);
};

} // namespace areko::remesher::extract
