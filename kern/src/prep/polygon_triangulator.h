#pragma once

#include "mesh/vec3.h"

#include <array>
#include <span>
#include <vector>

namespace areko::remesher::prep {

// Triangles of one polygon as triples of local corner numbers (0 .. n-1).
using LocalTriangles = std::vector<std::array<int, 3>>;

// Splits a polygon with any number of corners into n-2 triangles:
//  1. constrained Delaunay triangulation of the polygon projected onto its
//     plane (double); equal projected points are first moved apart a little;
//     taken only with exactly n-2 triangles and no new point,
//  2. otherwise the polygon tessellator on its own float projection; taken
//     only without new points and with exactly n-2 triangles,
//  3. otherwise a fan from corner 0.
class PolygonTriangulator {
public:
    [[nodiscard]] static LocalTriangles triangulate(std::span<const mesh::Vec3f> polygon);
    // Fan around corner `start`: (start, k, k+1) for every edge k -> k+1 not touching start.
    static void fan(int n, int start, LocalTriangles& out);
    // Float Newell normal (pairs i, i+1, the closing pair last), normalized
    // with the float normalization (shorter than 1e-28 stays as it is).
    [[nodiscard]] static mesh::Vec3f planeNormal(std::span<const mesh::Vec3f> polygon);
    // The same Newell sum in double, not normalized.
    [[nodiscard]] static mesh::Vec3d newellSum(std::span<const mesh::Vec3f> polygon);
    // Step 2 alone. With `rejectOnError` an error reported by the tessellator
    // rejects the result as well.
    [[nodiscard]] static bool tessellate(std::span<const mesh::Vec3f> polygon, LocalTriangles& out,
                                         bool rejectOnError = false);

private:
    [[nodiscard]] static bool delaunay(std::span<const mesh::Vec3f> polygon, LocalTriangles& out);
};

} // namespace areko::remesher::prep
