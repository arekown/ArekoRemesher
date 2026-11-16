#pragma once

#include "symmetry_cache.h"

#include "mesh/bit_array.h"
#include "mesh/mesh.h"

#include <span>

namespace areko::remesher::prep {

// Splits the long edges of a triangle mesh in passes. One pass:
//  - marks every edge longer than the target length (squared, float; with
//    protected borders no edge with a single corner),
//  - sharp corners (below 12 degrees) whose two edges are marked and whose
//    opposite edge is not get one point for both edges, unless the point's
//    valence, a fold of the edges (5 degrees or a locked edge) or the angles
//    in the neighbouring triangles forbid it,
//  - every other marked edge gets its midpoint (in edge order),
//  - the triangles are rebuilt in their order: split edges fan out, with two
//    split edges the shorter diagonal is taken (or, switched off, the fixed
//    diagonal from the split point after the apex),
//  - with mirror splitting every pass looks for the symmetry of the mesh as
//    it is; when found, the mirror edge of every long edge is split as well,
//    and without the shorter diagonal a triangle wholly on side 1 takes the
//    diagonal from the other split point.
// The vertex budget limits the points of the mesh; the pass that reaches it
// may add one point more. New points interpolate the point maps from their
// end points (sharp corner points by the weights of their position), face
// maps follow their faces, the sharp edges follow over point chains.
class EdgeSplitter {
public:
    // maxPasses -1 = unlimited, vertexBudget < 0 = unlimited.
    EdgeSplitter(int maxPasses, int vertexBudget, bool protectBorders) noexcept
        : maxPasses_(maxPasses), vertexBudget_(vertexBudget), protectBorders_(protectBorders)
    {
    }

    void setShortestDiagonal(bool on) noexcept { shortestDiagonal_ = on; }
    void setMirrorSplit(bool on) noexcept { mirrorSplit_ = on; }

    // Starts counting passes and points.
    void begin() noexcept;
    // One pass with the given target length; `locked` (may be null) marks
    // edges that count as folded for the sharp corner rule.
    void pass(mesh::Mesh& mesh, float targetLength, const mesh::BitArray* locked);
    // True while the passes added points and the maximum is not reached.
    [[nodiscard]] bool needsAnotherPass() const noexcept;

    [[nodiscard]] int passCount() const noexcept { return passCount_; }
    // Points of all passes; -1 before begin().
    [[nodiscard]] int insertedTotal() const noexcept { return insertedTotal_; }
    [[nodiscard]] int passInserted() const noexcept { return passInserted_; }
    // First point of the last pass and end of its sharp corner points.
    [[nodiscard]] int passFirstNew() const noexcept { return passFirstNew_; }
    [[nodiscard]] int passSharpEnd() const noexcept { return passSharpEnd_; }

    // Splits the listed edges at the given points (edge i at points[i]; an
    // edge listed twice keeps its first point) and rebuilds the triangles,
    // also when the list is empty. The map values of a new point come from the
    // end points (a, b) of its edge with the weights 1 - t[i] and t[i]. Two
    // split edges of a triangle take the shorter diagonal.
    static void splitAt(mesh::Mesh& mesh, std::span<const int> edges, std::span<const mesh::Vec3f> points,
                        std::span<const float> t);

    // Splits the listed edges at their midpoints (an edge listed twice keeps
    // its first point; with `mirror` pairs each edge is followed by its mirror
    // edge) and the triangles marked in `faces` (may be null) at their
    // centroid, then rebuilds the triangles. Without `shortestDiagonal` two
    // split edges take the fixed diagonal, or with `mirror` the one of the
    // side rule.
    static void splitMarked(mesh::Mesh& mesh, std::span<const int> edges, const mesh::BitArray* faces,
                            SymmetryCache* mirror, bool shortestDiagonal);

private:
    int maxPasses_;
    int vertexBudget_;
    bool protectBorders_;
    bool shortestDiagonal_ = true;
    bool mirrorSplit_ = false;
    int passCount_ = 0;
    int insertedTotal_ = -1;
    int passInserted_ = 0;
    int passFirstNew_ = -1;
    int passSharpEnd_ = -1;
};

} // namespace areko::remesher::prep
