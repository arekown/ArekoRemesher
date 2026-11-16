#pragma once

#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Reference length L of the preparation and the decimation length derived from it.
//  - frozen borders: the mean border edge length; with freeze marks the mean
//    chord between consecutive kept border points; above one million points
//    the normal rule with the quad count clamped to 3000 .. 12000 instead,
//  - otherwise the target edge length when > 0, else sqrt(area / quad count).
// With `areaCache` the area is taken from there while it is above 0 and
// stored there when it had to be computed.
class ReferenceLength {
public:
    [[nodiscard]] static float compute(const mesh::Mesh& mesh, float targetEdgeLength, int quadCount, bool freezeBorders,
                                       float* areaCache = nullptr);
    // Quad count clamped to 3000 .. 12000: below 8000 unchanged (at least
    // 3000), up to 14999 mapped linearly onto 8000 .. 12000 (float, truncated).
    [[nodiscard]] static int clampedQuadCount(int quadCount);
    // 1.2 L; with frozen borders at least sqrt(area / (2 * workResolution)).
    [[nodiscard]] static float decimationLength(const mesh::Mesh& mesh, float reference, bool freezeBorders,
                                                int workResolution);
    // Border length: robust mean of all border edges (two passes, the second
    // only with edges between 0.5 and 3 times the first mean), or with marks
    // the mean chord between consecutive kept points of every border loop
    // (FreezeMarks::keptPositions). Edges set in `locked` do not count, nor
    // does a chord that spans one. The mean is divided by `divisor` unless
    // that is 0; without any edge or chord the result is 0.
    [[nodiscard]] static float borderLength(const mesh::Mesh& mesh, const std::vector<mesh::FreezeMark>* marks,
                                            const mesh::BitArray* locked = nullptr, float divisor = 1.0f);

private:
    [[nodiscard]] static float normalRule(const mesh::Mesh& mesh, float targetEdgeLength, int quadCount,
                                          float* areaCache);
};

} // namespace areko::remesher::prep
