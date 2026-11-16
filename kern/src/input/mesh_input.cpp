#include "mesh_input.h"

#include "api/status.h"

namespace areko::remesher::input {

using api::StatusError;

void MeshInput::setMesh(std::span<const float> pointsXyz, std::span<const std::int32_t> faceSizes,
                        std::span<const std::int32_t> cornerPoints)
{
    *this = MeshInput{};
    if (faceSizes.empty()) {
        throw StatusError(ARM_ERROR_NO_MESH, "the mesh has no faces");
    }
    const auto points = static_cast<long long>(pointsXyz.size() / 3);
    long long corners = 0;
    for (const std::int32_t size : faceSizes) {
        if (size < 1) {
            throw StatusError(ARM_ERROR_INVALID_ARGUMENT, "a face has no corners");
        }
        corners += size;
    }
    if (corners != static_cast<long long>(cornerPoints.size())) {
        throw StatusError(ARM_ERROR_INVALID_ARGUMENT, "face sizes do not add up to the corner count");
    }
    for (const std::int32_t v : cornerPoints) {
        if (v < 0 || v >= points) {
            throw StatusError(ARM_ERROR_INVALID_ARGUMENT, "a corner names a point that does not exist");
        }
    }
    points_.assign(pointsXyz.begin(), pointsXyz.end());
    faceSizes_.assign(faceSizes.begin(), faceSizes.end());
    cornerPoints_.assign(cornerPoints.begin(), cornerPoints.end());
    hasMesh_ = true;
}

void MeshInput::requireMesh() const
{
    if (!hasMesh_) {
        throw StatusError(ARM_ERROR_NO_MESH, "no input mesh");
    }
}

void MeshInput::setFaceValues(arm_face_values kind, std::span<const std::int32_t> values)
{
    requireMesh();
    const auto index = static_cast<std::size_t>(kind);
    if (index >= faceValues_.size()) {
        throw StatusError(ARM_ERROR_INVALID_ARGUMENT, "unknown kind of face values");
    }
    faceValues_[index].emplace(values.begin(), values.end());
}

const std::optional<std::vector<std::int32_t>>& MeshInput::faceValues(arm_face_values kind) const
{
    return faceValues_.at(static_cast<std::size_t>(kind));
}

void MeshInput::setSharpEdges(std::span<const std::int32_t> pointPairs)
{
    requireMesh();
    for (const std::int32_t v : pointPairs) {
        if (v < 0 || v >= pointCount()) {
            throw StatusError(ARM_ERROR_INVALID_ARGUMENT, "a sharp edge names a point that does not exist");
        }
    }
    sharpEdges_.emplace(pointPairs.begin(), pointPairs.end());
}

void MeshInput::setCornerNormals(std::span<const float> normalsXyz, bool shareEqual)
{
    requireMesh();
    cornerNormals_.emplace(normalsXyz.begin(), normalsXyz.end());
    cornerNormalsShared_ = shareEqual;
}

void MeshInput::setCornerColors(std::span<const float> colorsRgba)
{
    requireMesh();
    cornerColors_.emplace(colorsRgba.begin(), colorsRgba.end());
}

void MeshInput::setPointDensity(std::span<const float> values)
{
    requireMesh();
    pointDensity_.emplace(values.begin(), values.end());
}

void MeshInput::addUvSet(std::string name, std::span<const float> uv)
{
    requireMesh();
    uvSets_.push_back({std::move(name), std::vector<float>(uv.begin(), uv.end())});
}

void MeshInput::addCurve(std::span<const float> pointsXyz, int kind, bool active, float width, float strength)
{
    requireMesh();
    if (kind < ARM_CURVE_BAND || kind > ARM_CURVE_OFF) {
        throw StatusError(ARM_ERROR_INVALID_ARGUMENT, "unknown curve kind");
    }
    if (pointsXyz.size() % 3 != 0) {
        throw StatusError(ARM_ERROR_INVALID_ARGUMENT, "curve points need three coordinates each");
    }
    curves_.push_back(Curve{std::vector<float>(pointsXyz.begin(), pointsXyz.end()), kind, active, width, strength});
}

} // namespace areko::remesher::input
