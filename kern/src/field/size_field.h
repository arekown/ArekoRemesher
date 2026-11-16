#pragma once

#include "mesh/bit_array.h"

#include <optional>
#include <vector>

namespace areko::remesher::field {

// Size field of one part: how long the sides of the quads are to be, per
// triangle along its two combed directions. Later steps only read it.
struct SizeField {
    // Target edge length of the part and the limits the sizes were kept
    // between while they followed the curvature.
    float edgeLength = 0.0f;
    float sizeMin = 0.0f;
    float sizeMax = 0.0f;
    // Per triangle the target size along the first and the second combed direction.
    std::vector<float> sizeU;
    std::vector<float> sizeV;
    // Number of quads these sizes amount to over the area of the part.
    int estimatedQuads = 0;

    // Only for a part with frozen borders.
    struct BorderData {
        // Per triangle the length of the path over triangle centres to the
        // nearest border triangle (FLT_MAX none), the triangle before it on
        // that path (itself at the border, -1 none) and the first corner of
        // the border edge the path starts at (-1 none).
        std::vector<float> distance;
        std::vector<int> previousFace;
        std::vector<int> borderCorner;
        // Border triangles whose size along the first resp. second direction
        // is given by the border and stays.
        mesh::BitArray fixedU;
        mesh::BitArray fixedV;
    };
    std::optional<BorderData> border;
};

} // namespace areko::remesher::field
