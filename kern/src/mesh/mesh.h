#pragma once

#include "adjacency.h"
#include "attribute_maps.h"
#include "connected_parts.h"
#include "normals.h"
#include "topology.h"
#include "uv_set.h"
#include "vec3.h"

#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace areko::remesher::mesh {

struct Bounds {
    Vec3f min;
    Vec3f max;
};

// Length of max - min (float) and the middle (min + max) * 0.5 per coordinate.
inline float diagonal(const Bounds& b) { return length(b.max - b.min); }
inline Vec3f center(const Bounds& b) { return (b.min + b.max) * 0.5f; }

// Points, faces, UV sets, maps and given normals of one mesh. Neighbourhood,
// parts, triangle geometry, edge angles and bounds are derived on first use
// and dropped when the data they depend on changes.
class Mesh {
public:
    Mesh() = default;
    // Copies data and state; derived caches are rebuilt on demand.
    Mesh(const Mesh& other) : Mesh(other, false) {}
    // With `keepCaches` the copy also takes over triangle geometry (including
    // replaced normals), edge angles and bounds; neighbourhood and parts are
    // always rebuilt on demand.
    Mesh(const Mesh& other, bool keepCaches);
    Mesh& operator=(const Mesh& other) = delete;
    Mesh(Mesh&& other) noexcept;
    Mesh& operator=(Mesh&& other) noexcept;

    // Exchanges the whole content of two meshes.
    void swap(Mesh& other) noexcept;
    // Takes over the content of `other`; this mesh keeps its replacement normal mode.
    void copyFrom(const Mesh& other);

    [[nodiscard]] int pointCount() const noexcept { return static_cast<int>(points_.size()); }
    [[nodiscard]] int faceCount() const noexcept { return topology_.faceCount(); }
    [[nodiscard]] int cornerCount() const noexcept { return topology_.cornerCount(); }

    [[nodiscard]] const std::vector<Vec3f>& points() const noexcept { return points_; }
    // Replaces all points (count may change; call topologyChanged() as well if faces refer to new points).
    void swapPoints(std::vector<Vec3f>& points);
    void appendPoints(std::span<const Vec3f> points);
    // Moves single points; drops every cache that depends on positions.
    void movePoints(std::span<const int> indices, std::span<const Vec3f> positions);

    [[nodiscard]] const Topology& topology() const noexcept { return topology_; }
    // Changes go through the returned reference and must be followed by topologyChanged().
    [[nodiscard]] Topology& editTopology() noexcept { return topology_; }
    void topologyChanged();
    void toPolygonMode();

    [[nodiscard]] AttributeMaps& maps() noexcept { return maps_; }
    [[nodiscard]] const AttributeMaps& maps() const noexcept { return maps_; }
    [[nodiscard]] std::vector<UvSet>& uvSets() noexcept { return uvSets_; }
    [[nodiscard]] const std::vector<UvSet>& uvSets() const noexcept { return uvSets_; }
    [[nodiscard]] std::optional<CustomNormals>& customNormals() noexcept { return customNormals_; }
    [[nodiscard]] const std::optional<CustomNormals>& customNormals() const noexcept { return customNormals_; }

    [[nodiscard]] const Adjacency& adjacency() const;
    [[nodiscard]] const ConnectedParts& connectedParts() const;
    [[nodiscard]] const TriangleGeometry& triangleGeometry() const;
    [[nodiscard]] const std::vector<float>& dihedralAngles() const;
    [[nodiscard]] const Bounds& bounds() const;

    // Mode of the replacement normals (0 = off). With apply = false only the
    // value changes; otherwise the replacement runs at once unless the same
    // mode was already applied.
    [[nodiscard]] int replacementNormalMode() const noexcept { return replacementMode_; }
    void setReplacementNormalMode(int mode, bool apply);

private:
    void geometryChanged();

    std::vector<Vec3f> points_;
    Topology topology_;
    AttributeMaps maps_;
    std::vector<UvSet> uvSets_;
    std::optional<CustomNormals> customNormals_;
    int replacementMode_ = 0;

    mutable std::unique_ptr<Adjacency> adjacency_;
    mutable std::unique_ptr<ConnectedParts> connectedParts_;
    mutable std::unique_ptr<TriangleGeometry> triangles_;
    mutable bool replacementApplied_ = false;
    mutable std::optional<std::vector<float>> dihedralAngles_;
    mutable std::optional<Bounds> bounds_;
};

} // namespace areko::remesher::mesh
