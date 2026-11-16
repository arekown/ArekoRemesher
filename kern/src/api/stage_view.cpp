#include "stage_view.h"

#include "mesh/sharp_edges.h"

namespace areko::remesher::api {

using namespace areko::remesher::mesh;

StageView::StageView(arm_stage stage, const Mesh& mesh, std::span<const prep::NetLine> lines,
                     std::span<const int> mirrorSides)
    : stage_(stage)
{
    points_.reserve(3 * mesh.points().size());
    for (const Vec3f& p : mesh.points()) {
        points_.insert(points_.end(), {p.x, p.y, p.z});
    }
    const Topology& topology = mesh.topology();
    faceStarts_.reserve(static_cast<std::size_t>(topology.faceCount()) + 1);
    faceStarts_.push_back(0);
    for (int f = 0; f < topology.faceCount(); ++f) {
        faceStarts_.push_back(topology.faceEnd(f));
    }
    corners_.assign(topology.cornerVertex().begin(), topology.cornerVertex().end());

    const AttributeMaps& maps = mesh.maps();
    if (const auto* m = maps.findAs<TypedMap<int>>(map_names::materials, ElementKind::Face)) {
        faceValues_[ARM_FACE_MATERIAL].emplace(m->values().begin(), m->values().end());
    }
    if (const auto* m = maps.findAs<TypedMap<int>>(map_names::faceGroups, ElementKind::Face)) {
        faceValues_[ARM_FACE_GROUP].emplace(m->values().begin(), m->values().end());
    }
    if (const auto* m = maps.findAs<TypedMap<std::uint32_t>>(map_names::smoothingGroups, ElementKind::Face)) {
        auto& out = faceValues_[ARM_FACE_SMOOTHING_GROUPS].emplace();
        for (const std::uint32_t v : m->values()) {
            out.push_back(static_cast<std::int32_t>(v));
        }
    }
    if (!mirrorSides.empty()) {
        faceValues_[ARM_FACE_MIRROR_SIDE].emplace(mirrorSides.begin(), mirrorSides.end());
    }
    const Adjacency& adjacency = mesh.adjacency();
    BitArray sharp = SharpEdges::edgeBits(mesh).value_or(BitArray{});
    sharp.resize(static_cast<std::size_t>(adjacency.edgeCount()));
    prep::NetLines::addEdges(mesh, lines, true, sharp);
    for (std::size_t e = 0; e < sharp.size(); ++e) {
        if (sharp.test(e)) {
            const auto ends = adjacency.edgeVertices(static_cast<int>(e));
            sharpEdgePairs_.insert(sharpEdgePairs_.end(), {ends[0], ends[1]});
        }
    }
    if (stage >= ARM_STAGE_FEATURE_LINES) {
        for (const prep::NetLine& line : lines) {
            linePoints_.insert(linePoints_.end(), line.vertices.begin(), line.vertices.end());
            lineStarts_.push_back(static_cast<std::int32_t>(linePoints_.size()));
        }
    }
    if (const auto& custom = mesh.customNormals()) {
        cornerNormal_.assign(custom->cornerNormal.begin(), custom->cornerNormal.end());
        for (const Vec3f& n : custom->normals) {
            normals_.insert(normals_.end(), {n.x, n.y, n.z});
        }
    }
    if (const auto* m = maps.findAs<TypedMap<Color4ub>>(map_names::colors, ElementKind::Point)) {
        for (const Color4ub& c : m->values()) {
            colors_.insert(colors_.end(), {c.r, c.g, c.b, c.a});
        }
    } else if (const auto* rgb = maps.findAs<TypedMap<Color3ub>>(map_names::colors, ElementKind::Point)) {
        for (const Color3ub& c : rgb->values()) {
            colors_.insert(colors_.end(), {c.r, c.g, c.b, std::uint8_t{255}});
        }
    }
    if (const auto* m = maps.findAs<TypedMap<float>>(map_names::density, ElementKind::Point)) {
        density_ = m->values();
    }
    for (const UvSet& set : mesh.uvSets()) {
        Uv uv{set.name, std::vector<std::int32_t>(set.cornerUv.begin(), set.cornerUv.end()), {}};
        for (const auto& c : set.coords) {
            uv.coords.insert(uv.coords.end(), {c[0], c[1]});
        }
        uvSets_.push_back(std::move(uv));
    }
}

void StageView::setParts(std::span<const int> faceParts, std::span<const int> partValues)
{
    if (!faceParts.empty()) {
        faceValues_[ARM_FACE_PART].emplace(faceParts.begin(), faceParts.end());
    }
    partValues_.assign(partValues.begin(), partValues.end());
}

void StageView::setLineKinds(std::span<const int> kinds)
{
    if (stage_ >= ARM_STAGE_FEATURE_LINES && kinds.size() + 1 == lineStarts_.size()) {
        lineKinds_.assign(kinds.begin(), kinds.end());
    }
}

void StageView::setFaceDirections(std::span<const float> directions)
{
    if (directions.size() == 3 * (faceStarts_.size() - 1)) {
        faceDirections_.assign(directions.begin(), directions.end());
    }
}

void StageView::setCrossField(std::span<const float> directions, std::span<const int> singularities)
{
    if (directions.size() == 3 * (faceStarts_.size() - 1)) {
        crossDirections_.assign(directions.begin(), directions.end());
    }
    singularities_.assign(singularities.begin(), singularities.end());
}

void StageView::setCutEdges(std::span<const int> pointPairs)
{
    cutEdgePairs_.assign(pointPairs.begin(), pointPairs.end());
}

void StageView::setSizeField(std::span<const float> faceSizes, std::span<const float> partSizes)
{
    if (faceSizes.size() == 2 * (faceStarts_.size() - 1)) {
        faceSizes_.assign(faceSizes.begin(), faceSizes.end());
        partSizes_.assign(partSizes.begin(), partSizes.end());
    }
}

void StageView::setGrid(std::span<const float> cornerCoordinates, std::span<const int> partSolved)
{
    if (cornerCoordinates.size() == 2 * corners_.size()) {
        gridCoordinates_.assign(cornerCoordinates.begin(), cornerCoordinates.end());
        gridParts_.assign(partSolved.begin(), partSolved.end());
    }
}

arm_mesh_view StageView::meshView() const noexcept
{
    arm_mesh_view view{};
    view.point_count = static_cast<int32_t>(points_.size() / 3);
    view.points_xyz = points_.data();
    view.face_count = static_cast<int32_t>(faceStarts_.size()) - 1;
    view.face_starts = faceStarts_.data();
    view.corner_count = static_cast<int32_t>(corners_.size());
    view.corner_points = corners_.data();
    return view;
}

const std::optional<std::vector<std::int32_t>>& StageView::faceValues(arm_face_values kind) const
{
    return faceValues_.at(static_cast<std::size_t>(kind));
}

} // namespace areko::remesher::api
