#pragma once

#include "cut_graph.h"
#include "smoothing_neighbours.h"

#include "mesh/adjacency.h"

#include <vector>

namespace areko::remesher::field {

// Smoothing of the two target sizes per triangle over the neighbours across
// inner edges. Across an edge with an odd number of quarter turns left between
// the combed crosses, the neighbour's two sizes change places. Across a curve
// edge only the size along the curve is taken over. Every filter reads the
// sizes as they were before its round and covers all triangles.
class SizeSmoothing {
public:
    SizeSmoothing(const mesh::Adjacency& adjacency, const CutGraph& cuts) noexcept
        : neighbours_(adjacency, cuts, SmoothingNeighbours::Limit::EveryCurve)
    {
    }

    // Mean of the own size and the sizes taken over from the neighbours.
    void mean(std::vector<float>& sizeU, std::vector<float>& sizeV) const;
    // `rounds` times the smallest of the own size and the sizes taken over.
    void minimum(int rounds, std::vector<float>& sizeU, std::vector<float>& sizeV) const;

private:
    SmoothingNeighbours neighbours_;
};

} // namespace areko::remesher::field
