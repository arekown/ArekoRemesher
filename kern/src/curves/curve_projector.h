#pragma once

#include "guide_curve.h"

#include "mesh/mesh.h"

#include <span>
#include <vector>

namespace areko::remesher::curves {

// Lays a located polyline onto a mesh (triangles or polygons) as pieces that
// run over its points, edges and faces. Face locations (locating, walks and
// paths over edges from a face) count three corners per face and hold for
// triangle meshes only; on a polygon mesh the polyline is located at points.
//  - Every point is located once more: a face location within the tolerance
//    of a corner point or of a side moves onto it (point before side).
//  - Two neighbouring locations are joined directly when they share a face,
//    an edge or a point; otherwise a straight walk in the plane of the faces
//    leads from one to the other, and the shortest path over edges is the
//    fallback (and the alternative when the walk went backwards: the shorter
//    one counts).
//  - Where neither reaches the next point the piece ends (kept with at least
//    two entries) and a new one starts.
//  - A piece that ends inside a face takes the last point; a closed curve
//    ends where it began.
// From a point the walk looks at every side of its faces that does not touch
// it and reaches a target point among the other corners of a face (across
// the face when the target is no neighbour of the point in it); from an edge
// at every other side of its faces.
// Tolerance: max(1e-6 diagonal, 1e-8) unless fixed; absolute tolerance 1e-5.
class CurveProjector {
public:
    struct Options {
        // > 0: fixed tolerance of the location, independent of the mesh size.
        float tolerance = 0.0f;
        // A straight walk is laid onto the chord between its two ends point by
        // point and given up when a point lies further from it than the mean
        // length of the edges at the ends that are points; a path over edges
        // has its points nearer than that mean moved onto the chord.
        bool chordCheck = false;
        // Every entry carries the number of the polyline point it reaches,
        // -1 for the entries between two of them.
        bool sourceCodes = false;
    };

    [[nodiscard]] static std::vector<CurveSegment> project(const mesh::Mesh& mesh, const LocatedPolyline& curve);
    [[nodiscard]] static std::vector<CurveSegment> project(const mesh::Mesh& mesh, const LocatedPolyline& curve,
                                                           const Options& options);
    // The same for several polylines on one unchanged mesh, in their order
    // (the face normals are computed once).
    [[nodiscard]] static std::vector<std::vector<CurveSegment>>
    projectAll(const mesh::Mesh& mesh, std::span<const LocatedPolyline> curves, const Options& options);

    // Shortest path over edges from point `from` to point `to` of the mesh:
    // the points after `from` up to and including `to`; empty without a path.
    [[nodiscard]] static std::vector<int> edgePath(const mesh::Mesh& mesh, int from, int to);
};

} // namespace areko::remesher::curves
