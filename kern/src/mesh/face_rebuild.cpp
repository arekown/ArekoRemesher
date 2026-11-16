#include "face_rebuild.h"

#include "mesh.h"
#include "modification.h"

namespace areko::remesher::mesh {

FaceRebuild::FaceRebuild(const Mesh& mesh, std::vector<int> cornerPoints)
    : topology_(mesh.topology()), cornerPoints_(std::move(cornerPoints))
{
}

int FaceRebuild::nextInFace(int f, int c) const noexcept
{
    return c + 1 < topology_.faceEnd(f) ? c + 1 : topology_.faceBegin(f);
}

int FaceRebuild::prevInFace(int f, int c) const noexcept
{
    return c > topology_.faceBegin(f) ? c - 1 : topology_.faceEnd(f) - 1;
}

void FaceRebuild::keep(int f, std::span<const int> corners)
{
    const std::size_t m = corners.size();
    if (m < 3) {
        reduced_ = true;
        return;
    }
    if (static_cast<int>(m) != topology_.faceSize(f)) {
        reduced_ = true;
    }
    for (std::size_t k = 0; k < m; ++k) {
        const int c = corners[k];
        const int following = corners[(k + 1) % m];
        int origin = -1;
        if (point(nextInFace(f, c)) == point(following)) {
            origin = c;
        } else if (const int before = prevInFace(f, following); point(before) == point(c)) {
            origin = before;
        }
        cornerVertex_.push_back(point(c));
        cornerOrigin_.push_back(c);
        halfEdgeOrigin_.push_back(origin);
    }
    faceStart_.push_back(static_cast<int>(cornerVertex_.size()));
    faceOrigin_.push_back(f);
}

void FaceRebuild::keepAll(int f)
{
    all_.clear();
    for (int c = topology_.faceBegin(f); c < topology_.faceEnd(f); ++c) {
        all_.push_back(c);
    }
    keep(f, all_);
}

void FaceRebuild::apply(Mesh& mesh, FaceMode mode)
{
    Modification change;
    change.faceStart = std::move(faceStart_);
    change.cornerVertex = std::move(cornerVertex_);
    change.faceMode = mode;
    change.points.state = 1;
    change.faces.newToOld = std::move(faceOrigin_);
    change.corners.newToOld = std::move(cornerOrigin_);
    change.halfEdges.newToOld = std::move(halfEdgeOrigin_);
    MeshModifier::apply(mesh, change);
}

} // namespace areko::remesher::mesh
