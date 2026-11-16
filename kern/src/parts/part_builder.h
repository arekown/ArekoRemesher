#pragma once

#include "part.h"

#include "mesh/mesh.h"
#include "prep/symmetry_cache.h"

#include <vector>

namespace areko::remesher::parts {

// How the maps of the meshes that PartBuilder::merged takes come together. By
// default (stage meshes) the maps of the first taken mesh grow with every
// later mesh: by the values of its same-named map, else by default values;
// edge maps do not grow. With `result` a taken mesh brings all its maps while
// the merged mesh has neither points nor faces, and edge maps grow by the
// edge count; with `unite` as well, a map that only a later mesh has is taken
// over, with default values for the elements before it.
struct MapRules {
    bool result = false;
    bool unite = false;
};

// Splits the prepared mesh into its connected pieces, gives each piece its
// share of the quad budget and decides which pieces are too small to compute.
class PartBuilder {
public:
    struct Split {
        std::vector<Part> parts;
        std::vector<int> faceParts;  // part of every face of the prepared mesh
        std::vector<int> facePlaces; // place of every face in its part
        bool sharedPoints = false;   // a point joins several pieces: no parts
    };

    // One part per component (seed: the smallest free face, breadth first),
    // in component order. A part copies the points of its component in
    // ascending order and its faces in face order, the point and face maps
    // with them; other maps, UV sets and normals stay behind. With symmetry on
    // and pairs of the cut, each part gets the pairs inside it (a partner in
    // another part becomes -1), the point sides and the locked plane.
    [[nodiscard]] static Split split(const mesh::Mesh& prepared, const prep::SymmetryCache* mirror, bool symmetryOn);

    // Quad budget per part (float): N = the budget, or with an edge length
    // the area over its square; one part takes trunc(N), several take
    // A_i * (N / A), at least 10, truncated. Every part keeps the edge length.
    static void assignBudgets(std::vector<Part>& parts, const mesh::Mesh& prepared, int quadBudget, float edgeLength);

    // Whether a part is too small to compute: two triangles on three points;
    // or by its share r of the prepared area and its point count (r < 2e-6:
    // fewer than 7 points, or up to 20 points with r < 1e-9; 2e-6 <= r < 1e-4
    // with at most 5 points: all edges steep (> 100 degrees) or on the border,
    // without border more than half of them steep).
    [[nodiscard]] static bool negligible(const Part& part, const mesh::Mesh& prepared);

    // The mesh of a part that a merge takes, or null when the part adds nothing.
    using MeshOf = const mesh::Mesh* (*)(const Part& part);
    // The raw mesh of a part (finished after the extraction), or null.
    [[nodiscard]] static const mesh::Mesh* rawMeshOf(const Part& part);

    // All parts one after another, each with the mesh `meshOf` gives it
    // (points, faces, the maps after `rules`); the faces stay triangles when
    // every taken mesh has triangles, else they are polygons. The part number
    // of every face goes to `faceParts`, the first point of every part to
    // `pointStarts` (one more entry at the end).
    [[nodiscard]] static mesh::Mesh merged(const std::vector<Part>& parts, MeshOf meshOf, std::vector<int>& faceParts,
                                           std::vector<int>& pointStarts, MapRules rules = {});
};

} // namespace areko::remesher::parts
