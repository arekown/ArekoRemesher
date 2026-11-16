#pragma once

#include "mesh/bit_array.h"
#include "mesh/mesh.h"

#include <vector>

namespace areko::remesher::prep {

// Repairs the topology of a mesh in this order:
//  1. split non-manifold points: every fan of a point after the first one
//     (fans join corners over edges with exactly two corners) gets a copy of
//     the point; first the points with the mask, then the others,
//  2. orient the faces part by part (breadth first from the first face of a
//     part, edges with more than two corners are skipped); a mesh that cannot
//     be oriented stays unchanged,
//  3. weld exactly equal points of each border loop onto the first of them in
//     x, y, z order (the first loop point stays itself), when all loop points
//     are equal or the welded loop stays manifold; afterwards degenerate
//     faces are cleaned,
//  4. split again while edges with more than two corners remain,
//  5. remove unused points.
// Every change keeps the maps of the mesh.
class MeshRepair {
public:
    // Returns removed unused points + new points from splitting + welded points.
    static int run(mesh::Mesh& mesh);
    // Drops points without corners (the rest keeps its order); returns their number.
    // `renumbering` (may be null) receives old -> new point numbers (-1 removed)
    // when points were removed.
    static int removeUnusedPoints(mesh::Mesh& mesh, std::vector<int>* renumbering = nullptr);
    // Step 1 alone, for the `candidates` only (null: all points; points beyond
    // its size do not count): the first pass leaves out the end point with
    // fewer corners of an isolated non-manifold edge, the second pass takes the
    // candidates the first one left out. Returns the source point of every
    // new point in order.
    [[nodiscard]] static std::vector<int> splitNonManifoldPoints(mesh::Mesh& mesh, const mesh::BitArray* candidates);
    // One pass of step 1 for the points set in `mask` (points beyond its size
    // do not count); appends the source point of every new point to `sources`.
    static void splitFans(mesh::Mesh& mesh, const mesh::BitArray& mask, std::vector<int>& sources);
    // Finishes a split: the mesh takes the corner list `cornerVertex` and one
    // new point per entry of `sources`, at the place of that point and with
    // its point map values; faces, corners and half edges stay as they are.
    static void splitPoints(mesh::Mesh& mesh, std::vector<int> cornerVertex, const std::vector<int>& sources);

private:
    static int splitNonManifold(mesh::Mesh& mesh);
    static bool orient(mesh::Mesh& mesh);
    static int weldLoopPoints(mesh::Mesh& mesh);
};

} // namespace areko::remesher::prep
