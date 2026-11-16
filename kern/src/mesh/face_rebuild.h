#pragma once

#include "topology.h"

#include <span>
#include <vector>

namespace areko::remesher::mesh {

class Mesh;

// New face list built from the faces of a mesh: every old face keeps some of
// its corners in their order or vanishes. The corners may name other points
// than before (one point per old corner). Applying the list keeps the points
// and their maps; faces, corners and half edges follow their origin. A new
// half edge takes the old half edge of the same face that joins the same two
// points: the one starting at its corner, else the one ending at the next
// kept corner, else none.
class FaceRebuild {
public:
    FaceRebuild(const Mesh& mesh, std::vector<int> cornerPoints);

    [[nodiscard]] int point(int corner) const noexcept { return cornerPoints_[static_cast<std::size_t>(corner)]; }

    // Keeps the listed old corners of face f (face order); with fewer than
    // three corners the face vanishes.
    void keep(int f, std::span<const int> corners);
    void keepAll(int f);

    // A face vanished or lost corners.
    [[nodiscard]] bool reduced() const noexcept { return reduced_; }

    // Sets the new faces; `mode` is the face mode afterwards when every face
    // fits it, polygons otherwise.
    void apply(Mesh& mesh, FaceMode mode);

private:
    [[nodiscard]] int nextInFace(int f, int c) const noexcept;
    [[nodiscard]] int prevInFace(int f, int c) const noexcept;

    const Topology& topology_;
    std::vector<int> cornerPoints_;
    std::vector<int> faceStart_{0};
    std::vector<int> cornerVertex_;
    std::vector<int> faceOrigin_;
    std::vector<int> cornerOrigin_;
    std::vector<int> halfEdgeOrigin_;
    std::vector<int> all_;
    bool reduced_ = false;
};

} // namespace areko::remesher::mesh
