#pragma once

#include "curvature.h"

#include "curves/guide_curve.h"
#include "mesh/mesh.h"

#include <vector>

namespace areko::remesher::field {

// What the guide curves of a part bring to the field steps: their width,
// the working polyline of bands, the tangents along it and the distance of
// every point to the curves.
class CurveInfluence {
public:
    // Active curves of kind band, edge loop or crease without a width (not
    // >= 0, NaN included) get min(0.25 * length, r_max) (minss), curves cut
    // into the mesh three quarters of it. The length is that of the curve's
    // points for bands and of the cut chain otherwise (float sum of the
    // pieces from 0).
    static void assignWidths(std::vector<curves::GuideCurve>& curves, float maxRadius);

    // Lays every active band onto `mesh` once: its source points become the
    // working polyline with the triangle under each (PolylineLocator without
    // snapping).
    static void prepareBands(std::vector<curves::GuideCurve>& curves, const mesh::Mesh& mesh);

    // Distance to the active curves (all but kind "off"), breadth first with a
    // queue that holds each point once; the bands are laid before:
    //  - seeds in curve order: the corners of the band triangles resp. the
    //    points of the cut chain, the first seed of a point wins (distance 0),
    //  - from a point v with distance d, seed point P: a neighbour w at
    //    a = |P[w] - P| gets c = 1.005 d when that is larger, else a, and takes
    //    it (with the seed) when c is below its distance and below twice the
    //    width of the seed's curve,
    //  - per triangle the mean of distance / width over its three corners
    //    ((d0/b0 + d1/b1) + d2/b2) / 3 when all three have a seed, else FLT_MAX.
    [[nodiscard]] static CurveDistance distance(const mesh::Mesh& mesh, const std::vector<curves::GuideCurve>& curves);

    // Unit tangents along a polyline: forward at the first point, backward at
    // the last, central in between; a single point gets zero.
    [[nodiscard]] static std::vector<mesh::Vec3f> tangents(const std::vector<mesh::Vec3f>& points);

private:
    // Sum of the piece lengths of a polyline.
    [[nodiscard]] static float length(const std::vector<mesh::Vec3f>& points);
};

} // namespace areko::remesher::field
