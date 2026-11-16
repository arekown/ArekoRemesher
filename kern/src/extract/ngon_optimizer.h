#pragma once

#include "mesh/mesh.h"

namespace areko::remesher::extract {

// How a part's polygon mesh is turned into quads.
struct NgonSettings {
    // Odd faces are paired first (with frozen border data), in up to 11 rounds.
    bool balanceParity = false;
    bool smoothNewPoints = false;
    int smoothingIterations = 2;
    float smoothingLambda = 0.5f;

    // With frozen border data: parity rounds, 4 smoothing steps of 0.7; else
    // with followed borders 2 smoothing steps of 0.5; else no smoothing.
    [[nodiscard]] static NgonSettings forPart(bool frozenBorderData, bool followBorders);
};

// Splits every face of five or more corners of a polygon mesh into quads
// (and single triangles where no quad fits), dissolves the doublets left
// between two quads once (edge dissolution without removing straight points;
// the doublet point stays unused) and smooths the new points. New faces
// follow the kept ones, new points follow the old ones.
class NgonOptimizer {
public:
    static void run(mesh::Mesh& mesh, const NgonSettings& settings);
};

} // namespace areko::remesher::extract
