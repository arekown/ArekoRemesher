#pragma once

#include <cstdint>
#include <vector>

namespace areko::remesher::mesh {

enum class FaceMode : int { Polygons = 0, Quads = 1, Triangles = 2 };

// Faces as one corner list. In polygon mode, face f owns corners
// faceStart[f] .. faceStart[f+1]-1; in quad / triangle mode it owns 4f..4f+3
// resp. 3f..3f+2 and faceStart is [0]. A corner is also the half edge from its
// point to the point of the next corner of the same face.
class Topology {
public:
    [[nodiscard]] FaceMode mode() const noexcept { return mode_; }
    [[nodiscard]] int faceCount() const noexcept;
    [[nodiscard]] int cornerCount() const noexcept { return static_cast<int>(cornerVertex_.size()); }
    [[nodiscard]] int faceBegin(int f) const noexcept;
    [[nodiscard]] int faceEnd(int f) const noexcept;
    [[nodiscard]] int faceSize(int f) const noexcept { return faceEnd(f) - faceBegin(f); }

    [[nodiscard]] const std::vector<int>& faceStart() const noexcept { return faceStart_; }
    [[nodiscard]] const std::vector<int>& cornerVertex() const noexcept { return cornerVertex_; }
    [[nodiscard]] int vertex(int c) const noexcept { return cornerVertex_[static_cast<std::size_t>(c)]; }

    // 1: all corner points >= 0, 2: some corner is -1 (deleted).
    [[nodiscard]] std::uint8_t indexState() const noexcept;

    // Swaps the given lists in (polygon mode); the arguments receive the old lists.
    void swapPolygons(std::vector<int>& faceStart, std::vector<int>& cornerVertex);
    // Swaps only the corner list and sets the mode; quad/triangle mode resets faceStart to [0].
    void swapCorners(std::vector<int>& cornerVertex, FaceMode mode);

    // Switches to quad or triangle mode when every face has 4 resp. 3 corners.
    bool trySetUniformMode(FaceMode mode);
    // Back to polygon mode (builds faceStart).
    void toPolygonMode();
    // No faces, polygon mode.
    void clear();

private:
    std::vector<int> faceStart_{0};
    std::vector<int> cornerVertex_;
    FaceMode mode_ = FaceMode::Polygons;
};

} // namespace areko::remesher::mesh
