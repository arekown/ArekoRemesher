#include "topology.h"

#include <algorithm>
#include <utility>

namespace areko::remesher::mesh {

int Topology::faceCount() const noexcept
{
    switch (mode_) {
    case FaceMode::Polygons:
        return static_cast<int>(faceStart_.size()) - 1;
    case FaceMode::Quads:
        return cornerCount() / 4;
    case FaceMode::Triangles:
        return cornerCount() / 3;
    }
    return 0;
}

int Topology::faceBegin(int f) const noexcept
{
    switch (mode_) {
    case FaceMode::Polygons:
        return faceStart_[static_cast<std::size_t>(f)];
    case FaceMode::Quads:
        return 4 * f;
    case FaceMode::Triangles:
        return 3 * f;
    }
    return 0;
}

int Topology::faceEnd(int f) const noexcept
{
    switch (mode_) {
    case FaceMode::Polygons:
        return faceStart_[static_cast<std::size_t>(f) + 1];
    case FaceMode::Quads:
        return 4 * f + 4;
    case FaceMode::Triangles:
        return 3 * f + 3;
    }
    return 0;
}

std::uint8_t Topology::indexState() const noexcept
{
    return std::any_of(cornerVertex_.begin(), cornerVertex_.end(), [](int v) { return v < 0; }) ? 2 : 1;
}

void Topology::swapPolygons(std::vector<int>& faceStart, std::vector<int>& cornerVertex)
{
    faceStart_.swap(faceStart);
    cornerVertex_.swap(cornerVertex);
    mode_ = FaceMode::Polygons;
}

void Topology::swapCorners(std::vector<int>& cornerVertex, FaceMode mode)
{
    cornerVertex_.swap(cornerVertex);
    mode_ = mode;
    if (mode_ != FaceMode::Polygons && faceStart_.size() >= 2) {
        faceStart_.assign(1, 0);
    }
}

bool Topology::trySetUniformMode(FaceMode mode)
{
    if (mode == FaceMode::Polygons) {
        toPolygonMode();
        return true;
    }
    const int size = mode == FaceMode::Quads ? 4 : 3;
    const int n = faceCount();
    for (int f = 0; f < n; ++f) {
        if (faceSize(f) != size) {
            return false;
        }
    }
    if (mode_ == FaceMode::Polygons && n > 0 && cornerCount() != n * size) {
        return false;
    }
    mode_ = mode;
    faceStart_.assign(1, 0);
    return true;
}

void Topology::toPolygonMode()
{
    if (mode_ == FaceMode::Polygons) {
        return;
    }
    const int n = faceCount();
    const int size = mode_ == FaceMode::Quads ? 4 : 3;
    faceStart_.resize(static_cast<std::size_t>(n) + 1);
    for (int f = 0; f <= n; ++f) {
        faceStart_[static_cast<std::size_t>(f)] = f * size;
    }
    mode_ = FaceMode::Polygons;
}

void Topology::clear()
{
    faceStart_.assign(1, 0);
    cornerVertex_.clear();
    mode_ = FaceMode::Polygons;
}

} // namespace areko::remesher::mesh
