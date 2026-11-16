#pragma once

#include "mesh/bit_array.h"
#include "mesh/mesh.h"

#include <span>
#include <vector>

namespace areko::remesher::prep {

// Laplace smoothing of point positions (Jacobi: every sweep reads the
// positions of the sweep before). The mesh itself stays unchanged; the caller
// takes over the returned positions.
//  - one sweep: p' = p (1 - w) + (sum of the edge neighbours) w / degree, all in
//    float and in edge order; points with weight 0 or without edges keep p,
//  - sweep 0 reads the mesh, later sweeps a copy of the previous result,
//  - fold guard (triangle meshes only; not for points with per point weight 0
//    or without edges, which keep p unchecked): a moved point that lies on the outer
//    side of the opposite edge of one of its triangles (seen with the given
//    triangle normal) is pushed 1 % of that edge inside; after two corrected
//    rounds it goes halfway back to its position before the sweep, after two
//    more it goes back completely and is smoothed again in up to
//    `maxRepairRuns` + 1 repair runs, which read the current result,
//  - when smoothing all points, the guard sees the neighbours of the result
//    under construction; when smoothing a list, it sees the source positions.
class Smoother {
public:
    struct Options {
        int iterations = 1;
        float uniformWeight = 1.0f;            // when no weight per point is given
        bool foldGuard = false;
        const mesh::BitArray* guardPoints = nullptr; // guard only these points (null: all)
        int maxRepairRuns = 4;
    };

    Smoother(const mesh::Mesh& mesh, const Options& options);

    // All points; `weights` has one value per point, otherwise the uniform weight applies.
    [[nodiscard]] std::vector<mesh::Vec3f> applyAll(std::span<const float> weights,
                                                    std::span<const mesh::Vec3f> triangleNormals) const;
    // Only the listed points (list order, all from the same source) with the uniform weight.
    [[nodiscard]] std::vector<mesh::Vec3f> applyTo(std::span<const int> points,
                                                   std::span<const mesh::Vec3f> triangleNormals) const;

private:
    struct Sweep;

    [[nodiscard]] std::vector<mesh::Vec3f> run(std::span<const float> weights, const std::span<const int>* list,
                                               std::span<const mesh::Vec3f> triangleNormals) const;

    const mesh::Mesh& mesh_;
    Options options_;
    bool guard_;
};

} // namespace areko::remesher::prep
