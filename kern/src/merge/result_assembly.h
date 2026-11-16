#pragma once

#include "mesh/mesh.h"
#include "parts/part.h"
#include "prep/preparation.h"

#include <vector>

namespace areko::remesher::merge {

// The one result mesh of a run, from the finished parts:
//  - with two parts or more, a part whose face group map is present but
//    empty gets a new group for all its faces: the smallest value from 0 that
//    no part has and that was not given before,
//  - the parts one after another (a part without raw mesh adds nothing), with
//    the maps of all of them (a map a part lacks gets default values for its
//    elements); the result is a polygon mesh,
//  - with the symmetry on and mirroring at the end, one cut per mirror axis
//    X, Y, Z at the plane position of the preparation: the negative side is
//    kept and mirrored, points on the plane are shared, sharp edges and face
//    maps follow, point maps only with the colours kept,
//  - the face groups become the host's materials, face groups and smoothing
//    groups again (when the preparation combined them),
//  - on request the points are split along the sharp edges, and the output is
//    cleaned (points between two edges dissolved, then non-manifold points
//    split; unused points stay),
//  - the point colours stay with the result,
//  - materials, face groups or smoothing groups without exactly one value
//    per face are dropped (the host takes them only with that length).
// Throws a StatusError (no result) when no face is left.
class ResultAssembly {
public:
    struct Options {
        bool keepColors = false;
        bool splitAtCreases = false;
        bool cleanOutput = false;
    };

    [[nodiscard]] static mesh::Mesh run(std::vector<parts::Part>& parts, const prep::Prepared& prepared,
                                        const Options& options);
};

} // namespace areko::remesher::merge
