#pragma once

#include "hole_closing.h"

#include "mesh/bit_array.h"
#include "mesh/mesh.h"
#include "mesh/vec3.h"

#include <vector>

namespace areko::remesher::parts {
struct Part;
} // namespace areko::remesher::parts

namespace areko::remesher::extract {

// The last relaxation of the quad mesh of a part (part.raw): three Jacobi
// steps that pull the points towards the edge lengths and directions of the
// size and cross field, each followed by a walk of every point over the
// triangle mesh of the part by its shift (SurfaceWalker), so the points stay
// on the surface and their sources follow.
//  - per edge a length range, a weight and a direction: from both ends, each
//    from the size field and the combed directions of its triangle (the
//    triangle of the point's source when that lies better along the edge),
//    else from the sizes around the source or the target edge length; along
//    a curve edge its own direction,
//  - per point the anisotropy of the curvature around its source, and its
//    mobility: 0 for a kept border point and for a point on a point of the
//    triangle mesh with three or more fixed edges, or two that bend by more
//    than about 45 degrees; fixed are the curve edges, the border edges and,
//    with the symmetry on a symmetric part, the edges on its mirror plane,
//  - lines in the quad mesh are the curve edges and the border edges; a point
//    on a line moves only along it, a point where lines end or cross stays.
class FinalRelaxation {
public:
    // Per edge of the quad mesh.
    struct EdgeTargets {
        std::vector<mesh::Vec3f> direction; // (0, 0, 0): no direction
        std::vector<float> lengthMin;
        std::vector<float> lengthMax;
    };

    // What one step reads besides the positions.
    struct StepInput {
        EdgeTargets edges;
        std::vector<float> anisotropy; // per point
        std::vector<float> mobility;   // per point
        mesh::BitArray lines;          // per edge of the quad mesh
    };

    // Relaxes part.raw (nothing without raw mesh or without the fields of the part).
    static void run(parts::Part& part, bool symmetry);

    // The edge targets of part.raw from the fields of the part (computed once, before the steps).
    [[nodiscard]] static EdgeTargets edgeTargets(parts::Part& part);

    // One Jacobi step over every point of `mesh` (positions of `mesh`); the
    // dihedral angles of its edges are read for points on lines. The values
    // of the last point with four edges carry over to the next one where a
    // part of the mix is weighted 0.
    [[nodiscard]] static std::vector<mesh::Vec3f> step(const mesh::Mesh& mesh, const StepInput& input);

    // Dihedral angle per edge of a polygon mesh: between the face normals
    // (Newell sums in double, normalized when longer than 1e-26) of an edge
    // with exactly two corners as in mesh::dihedralAngles (a NaN cosine gives
    // NaN), else 0.
    [[nodiscard]] static std::vector<float> dihedralAngles(const mesh::Mesh& mesh);
};

} // namespace areko::remesher::extract
