#pragma once

#include "adjacency.h"
#include "bit_array.h"
#include "csr.h"

#include <vector>

namespace areko::remesher::mesh {

// Connected parts over shared edges. Seeds are taken in ascending face order,
// each part grows breadth first (corners of a face in order, corners of an
// edge ascending), so part i is the i-th component in this order.
class ConnectedParts {
public:
    ConnectedParts(const Topology& topology, const Adjacency& adjacency);

    // Component number per face; edges set in `blocked` are not crossed (edges
    // beyond its size count as open). Components are numbered in the order of
    // their smallest face.
    [[nodiscard]] static std::vector<int> faceComponents(const Topology& topology, const Adjacency& adjacency,
                                                         const BitArray* blocked, int& componentCount);

    [[nodiscard]] int partCount() const noexcept { return partFaces_.rows(); }
    [[nodiscard]] int partOfFace(int f) const noexcept { return faceParts_[static_cast<std::size_t>(f)]; }
    [[nodiscard]] const std::vector<int>& faceParts() const noexcept { return faceParts_; }
    [[nodiscard]] const Csr& partFaces() const noexcept { return partFaces_; }
    // Part of a point: -1 none, -2 several parts.
    [[nodiscard]] int partOfPoint(int v) const noexcept { return pointParts_[static_cast<std::size_t>(v)]; }

private:
    std::vector<int> faceParts_;
    Csr partFaces_;
    std::vector<int> pointParts_;
};

} // namespace areko::remesher::mesh
