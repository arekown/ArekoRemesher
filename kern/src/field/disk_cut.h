#pragma once

#include "cut_copy.h"

#include "mesh/bit_array.h"
#include "run/progress.h"

#include <span>

namespace areko::remesher::field {

// Finds the cuts that are still missing to make a mesh that is already cut
// open into a disc.
//  - a search over the triangles of the copy builds a spanning forest: edges
//    at border points of the copy (unless the point is singular) and the
//    edges of the mirrored half cost nothing, every other edge costs one,
//  - curve edges are not crossed; when the search runs dry, one crossing of
//    the longest curve chain that still has an unreached triangle behind it
//    is opened,
//  - every inner edge that is not crossed by the forest is a candidate cut;
//    dangling ends are taken away, and two candidates of one triangle that
//    meet in a point of their own are replaced by the third side (up to
//    three rounds).
class DiskCut {
public:
    struct Input {
        const CutCopy& copy;
        const mesh::Adjacency& original;
        const mesh::BitArray& mirrorHalf;      // triangles; empty without symmetry
        std::span<const int> singularPoints;   // original points
        std::span<const int> curveEdges;       // original edges
    };

    // One bit per edge of the copy: the cuts to add.
    [[nodiscard]] static mesh::BitArray restCuts(const Input& input, const run::Progress& progress);
};

} // namespace areko::remesher::field
