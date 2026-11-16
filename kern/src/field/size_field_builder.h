#pragma once

#include "mesh/mesh.h"
#include "parts/part.h"
#include "run/progress.h"

#include <span>

namespace areko::remesher::field {

// Size field of a part: how long the sides of the quads are to be, per
// triangle along its two combed directions.
//  - the sizes start at the target edge length of the part and are cut down
//    where the surface bends, so that a full turn gets the wanted number of
//    points,
//  - mean and minimum filters spread the small sizes to their surroundings,
//    the side ratio is limited,
//  - with a strict quad count the field is computed again with adapted
//    values until the number of quads it amounts to fits the budget,
//  - a weighted least squares fit over all triangles makes the field smooth,
//  - frozen borders give the sizes beside them,
//  - at the end the density painted or given per point scales the sizes.
class SizeFieldBuilder {
public:
    struct Options {
        bool exactCount = false;           // keep the quad budget strictly
        float pointsPerCircle = 16.0f;     // points a full turn of the surface gets
        float densityRatio = 4.0f;         // how many times smaller than h a size may get; 1: no adaption
        float maxQuadRatio = 4.0f;         // largest ratio of the two sizes of a triangle
        float borderDensityInverse = 1.0f; // one over the density factor at frozen borders
        bool colorDensity = false;         // the point colours scale the sizes
        bool valueDensity = false;         // the density values of the points scale the sizes instead
    };

    // Writes the result to part.sizeField from the part's target edge length,
    // curvature and cut graph. A part without curvature, cross field or cut
    // graph gets none.
    static void compute(parts::Part& part, const Options& options, const run::Progress& progress);

    // Number of quads the sizes amount to: the mean of the area over the mean
    // quad area and of the sum of triangle area over quad area, truncated.
    [[nodiscard]] static int estimateQuads(const mesh::Mesh& mesh, std::span<const float> sizeU,
                                           std::span<const float> sizeV);
};

} // namespace areko::remesher::field
