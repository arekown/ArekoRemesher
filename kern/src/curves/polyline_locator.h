#pragma once

#include "closest_face.h"
#include "guide_curve.h"

#include "mesh/mesh.h"

#include <optional>
#include <span>

namespace areko::remesher::curves {

// Lays the points of a polyline onto a triangle mesh, point by point:
//  - the nearest face within a 40th of the diagonal of the mesh, the point
//    moves onto it; points without a face drop out,
//  - with snapping the nearest corner point of that face within a 3000th of
//    the diagonal instead (the point moves onto it),
//  - a polyline whose first and last point are equal counts as closed, a
//    single point included: its last point is not located, and the first
//    entry is repeated at the end when there is one.
// The face search is built with the first polyline that has points.
class PolylineLocator {
public:
    PolylineLocator(const mesh::Mesh& mesh, bool snapToCorners) : mesh_(mesh), snapToCorners_(snapToCorners) {}

    // Per located point its position, its location (a face, with snapping
    // also a point) and the face found for it.
    [[nodiscard]] LocatedPolyline locate(std::span<const mesh::Vec3f> points);

private:
    const mesh::Mesh& mesh_;
    bool snapToCorners_;
    std::optional<ClosestFaceSearch> search_;
};

} // namespace areko::remesher::curves
