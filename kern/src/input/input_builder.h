#pragma once

#include "areko_remesher.h"
#include "mesh/mesh.h"
#include "mesh_input.h"

namespace areko::remesher::input {

// Builds the input mesh of a run from the raw host data:
//  - points as given, faces in polygon mode, corners in the given order,
//  - UV sets: UV points = distinct (point, u, v), numbered ascending by
//    (64-bit word of v:u, point), -0 counts as +0,
//  - normals: distinct values, numbered ascending by the six 32-bit words of
//    the three doubles, -0 counts as +0,
//  - maps in this order: face groups, materials, smoothing groups, sharp edges
//    (edge bits, only with at least one pair), colors (only when density from
//    color or keeping colors is on), density,
//  - density: per point the factor its given value stands for (between 0.25
//    and 4; 1 for values below 0 and NaN),
//  - colors: the last corner of a point decides; channel = (int)(clamp01(v) * 255.0);
//    alpha 255 when density from color is on and colors are not kept;
//    points without corners stay (0, 0, 0, 255).
class InputBuilder {
public:
    [[nodiscard]] static mesh::Mesh build(const MeshInput& input, const arm_settings& settings);
};

} // namespace areko::remesher::input
