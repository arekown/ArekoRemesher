#pragma once

#include "hole_closing.h"

#include "mesh/bit_array.h"
#include "mesh/mesh.h"
#include "mesh/vec3.h"

#include <optional>
#include <span>
#include <vector>

namespace areko::remesher::extract {

// Moves a point over the surface of a triangle mesh instead of through
// space, so that it stays on the surface; its source (point, edge or
// triangle it lies on) is carried along.
//  - the shift is laid into the tangent plane at the source (normal of the
//    triangle; of a point or edge the float sum of the normals of its
//    triangles, normalized) and split into a length R and a direction, in
//    double,
//  - inside a triangle the point walks along the ray to the side it leaves
//    through (not the side it came in by, not a side parallel to the ray;
//    the later side wins a tie); when R
//    ends inside, it stops there; otherwise it crosses into the triangle
//    beyond, the direction laid into its plane (turned about the shared edge
//    when it nearly vanishes), and R shrinks by the distance (at least 1e-5),
//  - from a point or an edge it leaves through the fan of triangles around
//    it; on a fixed edge, and at a point with two fixed edges, it walks only
//    along them; a point with one or more than two fixed edges stays,
//  - a point that lies on the border of its triangle already takes that
//    corner or side as its new source without moving,
//  - it stops at a border, at an edge without two triangles, after 300 steps
//    or when R is used up.
// The fixed edges count only when there is one bit per edge of the mesh. A
// point without source (or with another source) does not move.
class SurfaceWalker {
public:
    SurfaceWalker(const mesh::Mesh& triangles, const mesh::BitArray& fixedEdges);

    // Walks `position` by `shift` from `source`; both are updated.
    void walk(PointSource& source, mesh::Vec3f& position, mesh::Vec3f shift);

private:
    // Crossing of the ray o + t d with the line a + s (b - a); none when both
    // are (nearly) parallel.
    struct Crossing {
        double t = 0.0; // along the ray
        double s = 0.0; // along the side, 0 .. 1 between a and b
    };
    [[nodiscard]] static std::optional<Crossing> ray(mesh::Vec3f a, mesh::Vec3f b, const mesh::Vec3d& o,
                                                     const mesh::Vec3d& d);
    [[nodiscard]] mesh::Vec3f normalOver(std::span<const int> corners) const;
    [[nodiscard]] bool fixed(int e) const { return useFixed_ && fixed_.test(static_cast<std::size_t>(e)); }
    // Direction d (float) turned from the plane with normal n1 into the plane with normal n2.
    [[nodiscard]] static mesh::Vec3f turned(mesh::Vec3f n1, mesh::Vec3f n2, mesh::Vec3f d);

    const mesh::Adjacency& adjacency_;
    const std::vector<mesh::Vec3f>& points_;
    const std::vector<mesh::Vec3f>& normals_;
    const mesh::BitArray& fixed_;
    bool useFixed_;
};

} // namespace areko::remesher::extract
