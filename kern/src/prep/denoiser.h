#pragma once

#include "feature_chains.h"
#include "prep_state.h"

#include "mesh/bit_array.h"
#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Feature preserving smoothing of a triangle mesh after a crease chain
// detection. Only points whose neighbourhood looks like a spike, a dent or
// noise (by angle defect and edge angles) move; flat areas, crease chains,
// fixed points, cube corners and the mesh border stay.
//  - crease chains (FeatureChains) on edge angles from normals in which
//    degenerate triangles take a neighbour's normal; those normals stay in the
//    mesh and serve as the old normals of the first round,
//  - rounds of: weight per point, weight smoothing, Laplace smoothing with
//    fold guard on the flat points, then a check for flipped triangles that
//    smooths their points once more or moves them back,
//  - denoise: two rounds, two plain weight smoothings, weights up to 0.6, the
//    end points of non-manifold edges stay as well,
//  - after the detection: four rounds, quiet and slightly bent points keep a
//    weight, a weight smoothing towards the largest neighbour weight, weights
//    up to 0.95.
class Denoiser {
public:
    enum class Mode { Denoise, AfterDetection };

    struct Options {
        Mode mode = Mode::Denoise;
        FeatureChains::Options chains;
        const mesh::BitArray* fixedPoints = nullptr; // used when it has one bit per point
    };

    // Detection and smoothing on a triangle mesh; returns the detected chains.
    static FeatureChains::Result run(mesh::Mesh& mesh, const Options& options);

    // The denoise step: fixed are the end points of face group borders and of
    // sharp edges. The caller decides whether to run it (strength above 0,
    // enough points).
    static void run(mesh::Mesh& mesh, const PrepState& state);
};

} // namespace areko::remesher::prep
