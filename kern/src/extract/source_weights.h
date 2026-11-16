#pragma once

#include "hole_closing.h"

#include "mesh/csr.h"
#include "mesh/mesh.h"
#include "mesh/vec3.h"

#include <span>
#include <vector>

namespace areko::remesher::extract {

// How every output point mixes from the points of the triangle mesh it lies
// on, by its source, all in float:
//  - on a point v: (v, 1),
//  - on an edge with end points A, B: t = ((P - A) . (B - A)) / |B - A|^2,
//    (A, 1 - t), (B, t),
//  - in a triangle with corners Q0, Q1, Q2: (u, v) from the normal equations
//    of P - Q0 over Q1 - Q0 and Q2 - Q0, (Q0, (1 - u) - v), (Q1, u), (Q2, v),
//  - without source (also beyond the end of `sources`): an empty row.
// Nothing is clamped, so a point outside its edge or triangle gets weights
// outside 0 .. 1. The triangle mesh must be in triangle mode.
struct SourceWeights {
    mesh::Csr rows;
    std::vector<float> weights; // parallel to rows.entries

    [[nodiscard]] static SourceWeights of(const mesh::Mesh& triangles, std::span<const mesh::Vec3f> points,
                                          std::span<const PointSource> sources);
};

} // namespace areko::remesher::extract
