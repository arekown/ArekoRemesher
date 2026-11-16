#pragma once

#include "hole_closing.h"

#include "mesh/mesh.h"

#include <span>

namespace areko::remesher::extract {

// Gives the output mesh the point colours of the triangle mesh: an RGB
// colour map stays RGB, an RGBA map RGBA; without one nothing happens. Every
// channel of an output point is the sum over its source weights (in row
// order, float) of weight times channel; a sum above 255 becomes 255, one
// below 0 becomes 0, the rest is truncated (a sum that is not a number gives
// 0). A point without source gets 0 in every channel.
class VertexColorTransfer {
public:
    static void run(mesh::Mesh& output, std::span<const PointSource> sources, const mesh::Mesh& triangles);
};

} // namespace areko::remesher::extract
