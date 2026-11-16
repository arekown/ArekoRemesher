#pragma once

#include "symmetry_cache.h"

#include "mesh/mesh.h"

#include <optional>

namespace areko::remesher::prep {

// Cuts a mesh at the plane x_axis = c and keeps the negative side; without
// `cutOnly` the mirrored copy of the kept half is added (faces i and h + i are
// mirror faces) and, with `buildCache`, the mirror pairs are returned.
//  - tolerance: from the bounding box diagonal d (d * 1e-6) unless given,
//  - with `splitEdgesFirst` the edges crossing the plane get a point on it
//    first (triangle meshes),
//  - plane: with `local` the bounding box middle, snapped to 0 when it lies
//    closer than the limit, otherwise `planePos` (`usePresetPlane`) or 0; the
//    limit is the tolerance for the edge split (the plane in
//    `Result::planePos`) and d * 1e-7 for the cut,
//  - points within the tolerance snap onto the plane; of an edge crossing the
//    plane the end chosen by the cross ratio moves onto it (straight, or along
//    its most distant neighbour across the plane),
//  - faces with more than two kept corners stay unless all of them lie on the
//    plane; point maps, face maps and sharp edges follow (origins as asked).
class SymmetryCut {
public:
    struct Job {
        int axis = 0;
        bool splitEdgesFirst = false;
        bool cutOnly = false;
        bool local = true;
        bool usePresetPlane = false;
        float planePos = 0.0f;
        bool buildCache = true;
        bool keepPointOrigin = true;
        bool keepFaceOrigin = true;
        bool keepHardEdges = true;
        float crossRatio = 0.2f;
        float tolerance = -1.0f; // < 0: from the bounding box
    };

    struct Result {
        float planePos = 0.0f; // plane of the edge split (0 without it)
        std::optional<SymmetryCache> cache;
    };

    [[nodiscard]] static Result run(mesh::Mesh& mesh, Job job);
};

} // namespace areko::remesher::prep
