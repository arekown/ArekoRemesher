#pragma once

#include "mesh/adjacency.h"
#include "mesh/bit_array.h"

#include <cstdint>
#include <span>

namespace areko::remesher::prep {

// Edges between faces that belong to different groups. The value of the face of
// the first corner of an edge is compared with the face of every further corner.
class FaceBorders {
public:
    // Edges with at least two corners whose faces carry different values; with
    // `withBoundary` also every edge with a single corner.
    [[nodiscard]] static mesh::BitArray valueBorders(const mesh::Adjacency& adjacency, std::span<const int> faceValues,
                                                     bool withBoundary);
    // Edges where the smoothing masks of the two faces share no bit (two faces
    // without any group therefore always form a border).
    [[nodiscard]] static mesh::BitArray smoothingBorders(const mesh::Adjacency& adjacency,
                                                         std::span<const std::uint32_t> faceMasks);
};

} // namespace areko::remesher::prep
