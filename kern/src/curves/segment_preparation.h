#pragma once

#include "guide_curve.h"

#include "mesh/mesh.h"

#include <span>
#include <vector>

namespace areko::remesher::curves {

// Tidies the pieces of the curves before they are cut into the mesh (one
// copy per piece, in order):
//  A. points on an edge closer than `snap` (parameter) to an end point move
//     onto it; the mesh point is remembered with the curve position (with
//     `snapToPoint` the piece takes the mesh point's position),
//  B. points between two corner points that lie on the edge of these two
//     points move onto its first end point,
//  C. a zig-zag over one edge (three points on it, the faces alternating)
//     becomes one point at the mean of the three,
//  D. in up to 5 rounds the neighbours of a corner point that lie on an edge
//     of that point move onto it,
// finally the remembered mesh points move to the curve positions.
class SegmentPreparation {
public:
    struct Options {
        float snap = 0.05f;
        bool snapToPoint = true;
    };

    [[nodiscard]] static std::vector<CurveSegment> prepare(mesh::Mesh& mesh, std::span<const CurveSegment* const> pieces,
                                                           const Options& options);

    // Whether a piece is closed: at least three entries, first and last on an
    // edge, either inside the same edge at the same parameter or on the same
    // corner point.
    [[nodiscard]] static bool isClosed(const mesh::Adjacency& adjacency, const CurveSegment& piece);

    // Mesh point of an entry on an edge at parameter 0 or 1, else -1.
    [[nodiscard]] static int cornerPoint(const mesh::Adjacency& adjacency, int edge, float t);
};

} // namespace areko::remesher::curves
