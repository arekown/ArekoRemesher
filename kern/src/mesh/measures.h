#pragma once

#include "adjacency.h"
#include "bit_array.h"
#include "vec3.h"

#include <vector>

namespace areko::remesher::mesh {

class Mesh;

// Geometric measures of triangle meshes. Every function computes from the
// current points; callers decide how long a result stays valid.

// Corner angles per corner (radians, float). For triangle (P0, P1, P2):
// u = P2 - P0, w = P1 - P0, z = P1 - P2 in float; a0 = acos(u.w / (|u| |w|)) and
// a1 = acos(w.z / (|w| |z|)) with the dot product in float, lengths and
// quotient in double (denominator 0 gives 0, quotient <= -1 gives pi, >= 1
// gives 0); a2 = (pi - a0) - a1 in double from the unrounded a0, a1.
// Empty for meshes that are not in triangle mode.
[[nodiscard]] std::vector<float> cornerAngles(const Mesh& mesh);

// Angle between two vectors: zero when |u| |v| is 0 (lengths from float square
// sums, root and product in double), else acos(u.v / (|u| |v|)) with the dot
// product in float and the quotient in double (pi as float for <= -1, 0 for >= 1).
[[nodiscard]] float vectorAngle(Vec3f u, Vec3f v);

// Angle between the lines along u and v (0 .. pi / 2): zero when |u| |v| is 0
// (lengths and product as in vectorAngle), else acos of |u.v| / (|u| |v|) with
// the dot product in float and the quotient in double, 0 for >= 1.
[[nodiscard]] float lineAngle(Vec3f u, Vec3f v);

// Angle from u to v in 0 .. 2 pi, in double: the angle a between them as in
// vectorAngle (acos in double), taken as 2 pi - a unless the side
// crossed . normal (float) is positive; a = 0 stays 0, quotient <= -1 gives pi.
[[nodiscard]] double signedVectorAngleAsDouble(Vec3f u, Vec3f v, Vec3f crossed, Vec3f normal);

// The same angle rounded to float.
[[nodiscard]] float signedVectorAngle(Vec3f u, Vec3f v, Vec3f crossed, Vec3f normal);

// Angle at one corner of a triangle mesh between the edges to its previous and
// next corner (vectorAngle).
[[nodiscard]] float cornerAngle(const Mesh& mesh, int corner);

// Angle defect per point: 2 pi minus the sum of its (float) corner angles in
// corner order, both in double, rounded to float; points without corners keep
// 2 pi (float), border points get 0.
[[nodiscard]] std::vector<float> angleDefects(const Adjacency& adjacency, const std::vector<float>& cornerAngles);

// Angle defect per point of a triangle mesh without stored corner angles: per
// corner of the point (corner order) the angle between the edges to the next
// and the previous corner in double (as vectorAngle, without rounding to
// float), summed in double, 2 pi minus the sum rounded to float; points without
// corners keep 2 pi (float), border points get 0. It differs from angleDefects
// in the last bits.
[[nodiscard]] std::vector<float> directAngleDefects(const Mesh& mesh);

// Length of every edge (edgeLength of its two end points: in double, kept as float).
[[nodiscard]] std::vector<float> edgeLengths(const Mesh& mesh);

// Sets both end points of every set edge in `points`; the mask grows to
// `pointCount` first and keeps its bits.
void addEdgeEndPoints(const Adjacency& adjacency, const BitArray& edges, int pointCount, BitArray& points);

} // namespace areko::remesher::mesh
