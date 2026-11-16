#pragma once

#include "guide_curve.h"

#include "mesh/mesh.h"

#include <span>
#include <vector>

namespace areko::remesher::curves {

// Cuts curve pieces into a mesh (triangle, quad or polygon mesh), so that a
// chain of mesh edges follows every piece:
//  A. every entry becomes a mesh point: corner points stay, entries inside an
//     edge or a face become new points interpolated from the edge ends or the
//     face corners,
//  B./C. every face crossed by a piece gets an entry per crossing, given as
//     positions on its border (corner k = k, side k at parameter t = k + t)
//     or inside it,
//  D. new points at the same place on a border or inside are merged,
//  E. faces with entries are rebuilt: with an inside point by a constrained
//     Delaunay triangulation of the triangle, otherwise by splitting the
//     polygon along its chords,
//  F. on a triangle mesh the new polygons are triangulated,
//  G. the change is applied only when something changed,
//  H. two chain points without an edge between them are joined over a new
//     common neighbour or over the shortest edge path,
//  I. optionally unused points are removed.
// Every new point gets exactly one row of interpolation sources.
class SegmentCutter {
public:
    struct Options {
        bool compactPoints = true;   // remove unused points afterwards
        bool moveCornerPoints = false; // corner points hit by a piece move to its position
    };

    struct Result {
        std::vector<std::vector<int>> chains;  // mesh points per piece
        std::vector<std::vector<int>> numbers; // entry of the piece per chain point
        bool changed = false;                  // the mesh changed
    };

    [[nodiscard]] static Result cut(mesh::Mesh& mesh, std::span<const CurveSegment> pieces, const Options& options);
};

} // namespace areko::remesher::curves
