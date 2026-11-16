#include "mesh.h"

#include <algorithm>
#include <utility>

namespace areko::remesher::mesh {

Mesh::Mesh(const Mesh& other, bool keepCaches)
    : points_(other.points_),
      topology_(other.topology_),
      maps_(other.maps_),
      uvSets_(other.uvSets_),
      customNormals_(other.customNormals_),
      replacementMode_(other.replacementMode_)
{
    if (keepCaches) {
        if (other.triangles_) {
            triangles_ = std::make_unique<TriangleGeometry>(*other.triangles_);
        }
        replacementApplied_ = other.replacementApplied_;
        dihedralAngles_ = other.dihedralAngles_;
        bounds_ = other.bounds_;
    }
}

Mesh::Mesh(Mesh&& other) noexcept
{
    swap(other);
}

Mesh& Mesh::operator=(Mesh&& other) noexcept
{
    if (this != &other) {
        Mesh empty;
        swap(empty);
        swap(other);
    }
    return *this;
}

void Mesh::swap(Mesh& other) noexcept
{
    using std::swap;
    swap(points_, other.points_);
    swap(topology_, other.topology_);
    swap(maps_, other.maps_);
    swap(uvSets_, other.uvSets_);
    swap(customNormals_, other.customNormals_);
    swap(replacementMode_, other.replacementMode_);
    swap(adjacency_, other.adjacency_);
    swap(connectedParts_, other.connectedParts_);
    swap(triangles_, other.triangles_);
    swap(replacementApplied_, other.replacementApplied_);
    swap(dihedralAngles_, other.dihedralAngles_);
    swap(bounds_, other.bounds_);
    // The neighbourhood refers to the topology object, whose content just moved.
    if (adjacency_) {
        adjacency_->attach(topology_);
    }
    if (other.adjacency_) {
        other.adjacency_->attach(other.topology_);
    }
}

void Mesh::copyFrom(const Mesh& other)
{
    if (this == &other) {
        return;
    }
    const int mode = replacementMode_;
    Mesh copy(other);
    swap(copy);
    replacementMode_ = mode;
}

void Mesh::swapPoints(std::vector<Vec3f>& points)
{
    points_.swap(points);
    geometryChanged();
}

void Mesh::appendPoints(std::span<const Vec3f> points)
{
    points_.insert(points_.end(), points.begin(), points.end());
    geometryChanged();
}

void Mesh::movePoints(std::span<const int> indices, std::span<const Vec3f> positions)
{
    const std::size_t n = std::min(indices.size(), positions.size());
    for (std::size_t i = 0; i < n; ++i) {
        points_[static_cast<std::size_t>(indices[i])] = positions[i];
    }
    geometryChanged();
}

void Mesh::geometryChanged()
{
    triangles_.reset();
    replacementApplied_ = false;
    dihedralAngles_.reset();
    bounds_.reset();
}

void Mesh::topologyChanged()
{
    adjacency_.reset();
    connectedParts_.reset();
    geometryChanged();
}

void Mesh::toPolygonMode()
{
    if (topology_.mode() != FaceMode::Polygons) {
        topology_.toPolygonMode();
        topologyChanged();
    }
}

const Adjacency& Mesh::adjacency() const
{
    if (!adjacency_ || adjacency_->pointCount() != pointCount()) {
        adjacency_ = std::make_unique<Adjacency>(topology_, pointCount());
        connectedParts_.reset(); // the parts are derived from the neighbourhood
    }
    return *adjacency_;
}

const ConnectedParts& Mesh::connectedParts() const
{
    const Adjacency& neighbourhood = adjacency();
    if (!connectedParts_) {
        connectedParts_ = std::make_unique<ConnectedParts>(topology_, neighbourhood);
    }
    return *connectedParts_;
}

const TriangleGeometry& Mesh::triangleGeometry() const
{
    if (!triangles_) {
        triangles_ = std::make_unique<TriangleGeometry>(*this);
        replacementApplied_ = false;
        if (topology_.mode() == FaceMode::Triangles && replacementMode_ != 0) {
            triangles_->replaceDegenerateNormals(*this, adjacency(), replacementMode_);
            replacementApplied_ = true;
        }
    }
    return *triangles_;
}

void Mesh::setReplacementNormalMode(int mode, bool apply)
{
    const int previous = replacementMode_;
    replacementMode_ = mode;
    if (!apply || (previous == mode && triangles_ && replacementApplied_) || mode == 0) {
        return;
    }
    // Making the triangle geometry valid may already replace once; the explicit
    // request replaces again on the current normals.
    (void)triangleGeometry();
    triangles_->replaceDegenerateNormals(*this, adjacency(), mode);
}

const std::vector<float>& Mesh::dihedralAngles() const
{
    if (!dihedralAngles_) {
        dihedralAngles_ = mesh::dihedralAngles(adjacency(), triangleGeometry().normals());
    }
    return *dihedralAngles_;
}

const Bounds& Mesh::bounds() const
{
    if (!bounds_) {
        Bounds b;
        if (!points_.empty()) {
            b.min = b.max = points_.front();
            for (const Vec3f& p : points_) {
                b.min = {std::min(b.min.x, p.x), std::min(b.min.y, p.y), std::min(b.min.z, p.z)};
                b.max = {std::max(b.max.x, p.x), std::max(b.max.y, p.y), std::max(b.max.z, p.z)};
            }
        }
        bounds_ = b;
    }
    return *bounds_;
}

} // namespace areko::remesher::mesh
