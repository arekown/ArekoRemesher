#pragma once

#include "mesh/mesh.h"

#include <optional>
#include <vector>

namespace areko::remesher::field {

// Density factor per point of a mesh: 1 leaves the size of the quads as it
// is, 4 makes their sides four times shorter, 0.25 four times longer.
class DensityFactors {
public:
    // From the point colours: white 1, red up to 4, cyan down to 0.25
    // (factor = (1 + 3 (1 - (g + b) / 510)) / (1 + 3 (1 - r / 255)), float).
    // None without colours or when their number is not the point count.
    [[nodiscard]] static std::optional<std::vector<float>> fromColors(const mesh::Mesh& mesh);
    // The factor a density value stands for: the value itself between 0.25
    // and 4, the nearer limit outside, 1 for values below 0 and NaN. A factor
    // stays what it is.
    [[nodiscard]] static float fromDensity(float density);
    // From the density map of the points. The input stores the factors of
    // the given values there and the preparation mixes them like the colours;
    // a mixed value that is no number counts as 1. None without the map or
    // when its size is not the point count.
    [[nodiscard]] static std::optional<std::vector<float>> fromValues(const mesh::Mesh& mesh);

    // Divides both sizes of every triangle by the mean factor of its corners.
    static void apply(const mesh::Mesh& mesh, const std::vector<float>& factors, std::vector<float>& sizeU,
                      std::vector<float>& sizeV);
};

} // namespace areko::remesher::field
