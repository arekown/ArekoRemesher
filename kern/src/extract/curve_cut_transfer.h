#pragma once

#include "curve_marks.h"
#include "raw_mesh.h"

#include "mesh/adjacency.h"
#include "mesh/bit_array.h"
#include "mesh/mesh.h"

#include <span>
#include <vector>

namespace areko::remesher::extract {

// The output points on one crease chain in chain order: their numbers and
// their place along the chain (edge place plus the share of the edge).
struct CurvePolyline {
    std::vector<int> points;
    std::vector<float> keys;
};

// Cuts the crease chains into the quad mesh of a part:
//  - every polyline (its output points) is laid onto the quad mesh with a
//    fixed location tolerance, chord checks and the number of the polyline
//    point per reached entry,
//  - the entries between two reached points move onto the crease chain of
//    the triangle mesh by their share of the arc length and take the chain
//    edge (or end point) there as their source,
//  - the pieces are tidied and cut into the mesh; the corner points they hit
//    move onto them, unused points stay,
//  - the point data grow with the new points (no source, weight 0, on the
//    curve, no triangle, no neighbours); a new chain point whose entry has a
//    source takes it,
//  - the edges between consecutive chain points become the curve edges.
class CurveCutTransfer {
public:
    static void run(const mesh::Mesh& triangles, const CurveMarks& marks, std::span<const CurvePolyline> polylines,
                    RawMesh& raw);

    // Sets the edge between every two consecutive points of `chain` (where
    // there is one) in `edges`.
    static void markEdges(const mesh::Adjacency& adjacency, std::span<const int> chain, mesh::BitArray& edges);
};

} // namespace areko::remesher::extract
