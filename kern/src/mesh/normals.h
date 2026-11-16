#pragma once

#include "adjacency.h"
#include "vec3.h"

#include <cstdint>
#include <span>
#include <vector>

namespace areko::remesher::mesh {

class Mesh;

// Normals given with the input: every corner names one normal.
struct CustomNormals {
    std::vector<int> cornerNormal;
    std::vector<Vec3f> normals;
    std::uint8_t state = 1; // 2 when a corner has no normal (-1)

    [[nodiscard]] bool empty() const noexcept { return cornerNormal.empty(); }
};

// Unit normal (float, computed in double) and doubled area of the triangle (a, b, c).
void floatTriangleNormal(Vec3f a, Vec3f b, Vec3f c, Vec3f& normal, float& doubleArea);

// Centre of every triangle of a triangle mesh: per coordinate the three
// corner points summed in corner order and divided by 3 in double, then
// rounded to float.
[[nodiscard]] std::vector<Vec3f> triangleCentres(const Mesh& mesh);

// Per-triangle geometry: unit normals and doubled areas (cross product in
// double), the total area as float sum in triangle order times 0.5.
// Triangles: in triangle mode face t; other modes use a fan per face.
class TriangleGeometry {
public:
    explicit TriangleGeometry(const Mesh& mesh);

    [[nodiscard]] const std::vector<Vec3f>& normals() const noexcept { return normals_; }
    [[nodiscard]] std::vector<Vec3f>& normals() noexcept { return normals_; }
    [[nodiscard]] const std::vector<float>& doubleAreas() const noexcept { return doubleAreas_; }
    [[nodiscard]] float totalArea() const noexcept { return totalArea_; }
    // Limit of the doubled area of tiny triangles: float(double(total area) *
    // 1e-8), at least 1e-12 as float.
    [[nodiscard]] float smallAreaLimit() const noexcept;

    // Replacement normals for degenerate triangles (triangle meshes only, mode 2):
    // a triangle whose doubled area is not above a limit takes the normal of a
    // valid neighbour, in up to 20 rounds, replaced ones count as valid at once.
    void replaceDegenerateNormals(const Mesh& mesh, const Adjacency& adjacency, int mode);

    // Brings the values up to date after points moved (triangle meshes): every
    // triangle at one of `points` gets its normal and doubled area again,
    // without replacement; all other triangles keep their values. The total
    // area is summed anew.
    void updateAt(const Mesh& mesh, std::span<const int> points);

private:
    std::vector<Vec3f> normals_;
    std::vector<float> doubleAreas_;
    float totalArea_ = 0.0f;
};

// One normal per corner group around each point. Given normals that name a
// normal for every corner are taken as they are. Otherwise, on triangle
// meshes, the fan of every point (orderedVertexFan) is cut into groups: the
// next triangle joins the current group unless the cosine between its normal
// and the previous triangle's normal is below a limit of about 50 degrees or
// (with a sharp edge map per half edge) the half edge of its corner is sharp.
// Last and first triangle of the fan are compared the same way (with the
// sharp bit of the first corner); when they join, the normal of the last
// group is dropped and its corners take the group of the first corner. A
// group's normal is the sum of triangle normal times doubled area in fan
// order, normalized at the end; corners outside the fan share one (1, 0, 0)
// per point.
class SplitNormals {
public:
    explicit SplitNormals(const Mesh& mesh);

    [[nodiscard]] Vec3f atCorner(int corner) const;
    // Normalized sum, over the corners of edge e in order, of the normals at
    // the corner and at its next corner; zero for edges with more than two
    // corners.
    [[nodiscard]] Vec3f edgeNormal(const Adjacency& adjacency, int e) const;

private:
    void build(const Mesh& mesh);

    std::vector<int> cornerNormal_;
    std::vector<Vec3f> normals_;
};

// Angle between the normals of the two faces at each edge (one normal per
// face, any face mode): the float dot product, its angle as in angleOfCosine
// (NaN stays NaN); boundary and non-manifold edges get 0.
[[nodiscard]] std::vector<float> dihedralAngles(const Adjacency& adjacency, const std::vector<Vec3f>& faceNormals);

} // namespace areko::remesher::mesh
