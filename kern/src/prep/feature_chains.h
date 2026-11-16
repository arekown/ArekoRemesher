#pragma once

#include "mesh/bit_array.h"
#include "mesh/csr.h"
#include "mesh/mesh.h"

#include <memory>
#include <span>
#include <vector>

namespace areko::remesher::prep {

// Chains of crease edges on a triangle mesh and their class.
//  - candidates: edges outside the protected set whose edge angle (the mesh's
//    cached edge angles) is not below the crease angle,
//  - in rounds: the candidates are split into chains (EdgeChains), every chain
//    gets 19 measures, a first pass sets a preliminary class, a second pass
//    decides; rejected chains leave the candidates, some of them are kept
//    aside as rescued chains. Points where a chain breaks only lower its
//    smoothness, they do not split it. A round without rejection ends; after
//    ten rounds the next chains all get class 0,
//  - result: the chains of the last round, then the rescued ones (class 3).
// Classes: -1 rejected, 0 undecided, 1 and 2 crease lines, 3 soft line.
class FeatureChains {
public:
    // Reference lengths from area A and quad count N, in float:
    // reference = sqrt(A / D) with D = N raised to 1000 and limited to 10000,
    // quad = sqrt(A / max(N, 1)). A negative area gives NaN; every comparison
    // with them is written so that NaN takes a fixed branch.
    struct Lengths {
        float reference = 0.0f;
        float quad = 0.0f;
    };
    [[nodiscard]] static Lengths referenceLengths(float area, int quadCount);

    struct Options {
        float creaseAngle = 0.5235988f;                 // radians
        Lengths lengths;
        const mesh::BitArray* protectedEdges = nullptr; // never candidates (null: none)
    };

    struct Result {
        mesh::Csr points;         // chain -> points (a closed chain repeats its first point)
        mesh::Csr edges;          // chain -> edges
        std::vector<int> classes; // one per chain
    };

    [[nodiscard]] static Result classify(const mesh::Mesh& mesh, const Options& options);

    // Judges lines on one mesh by the same chain measures without candidates
    // or protected edges, measured once for all lines. The mesh must outlive
    // the check and keep its geometry.
    class CreaseLineCheck {
    public:
        explicit CreaseLineCheck(const mesh::Mesh& mesh);
        ~CreaseLineCheck();
        CreaseLineCheck(const CreaseLineCheck&) = delete;
        CreaseLineCheck& operator=(const CreaseLineCheck&) = delete;

        // Whether a line through the given points stands for a crease: its
        // smallest edge angle is above 45 degrees, or its largest and its
        // length weighted edge angle are not below 30 degrees while its bend
        // stays below 8 and its edge angle jumps below 10 degrees (limits in
        // double). A line with fewer than two points or two consecutive
        // points without an edge is no crease.
        [[nodiscard]] bool operator()(std::span<const int> points) const;

    private:
        class Measure;
        std::unique_ptr<const Measure> measure_;
    };
};

} // namespace areko::remesher::prep
