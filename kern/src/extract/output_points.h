#pragma once

#include "mesh/bit_array.h"

#include <array>
#include <cstddef>
#include <vector>

namespace areko::remesher::extract {

// Where an output point comes from: type 1 a point, 2 an edge, 4 a triangle
// of the triangle mesh (id = its number), 0 none.
struct PointSource {
    int type = 0;
    int id = -1;
    friend bool operator==(const PointSource&, const PointSource&) = default;
};

// Data per point of the raw mesh. The bit sets and the triangle and
// neighbour lists may be shorter than the point count (missing entries count
// as not set); source and weight cover every point.
struct OutputPoints {
    mesh::BitArray onBorderEdge;  // on a border edge of the triangle mesh
    mesh::BitArray reachedBorder; // on the border or walked into it
    std::vector<PointSource> source;
    std::vector<float> weight;
    std::vector<int> triangle;
    std::array<std::vector<int>, 4> neighbours;
    mesh::BitArray onCurve; // on a crease curve (set when the features are transferred)

    // Brings every list to `count` points (longer or shorter); new points get
    // no bits, `source` and `weight`, no triangle and no neighbours (-1).
    void resize(std::size_t count, PointSource newSource, float newWeight);
};

} // namespace areko::remesher::extract
