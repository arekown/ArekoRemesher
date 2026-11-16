#pragma once

#include "areko_remesher.h"

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace areko::remesher::input {

// Raw input of one run, as handed over by the host. Checked on arrival;
// turned into a mesh by InputBuilder.
class MeshInput {
public:
    // Starts a new input; every optional part is dropped.
    void setMesh(std::span<const float> pointsXyz, std::span<const std::int32_t> faceSizes,
                 std::span<const std::int32_t> cornerPoints);
    void setFaceValues(arm_face_values kind, std::span<const std::int32_t> values);
    void setSharpEdges(std::span<const std::int32_t> pointPairs);
    // shareEqual: equal normals get one index; otherwise every corner keeps its own.
    void setCornerNormals(std::span<const float> normalsXyz, bool shareEqual);
    void setCornerColors(std::span<const float> colorsRgba);
    void setPointDensity(std::span<const float> values);
    void addUvSet(std::string name, std::span<const float> uv);
    // One guide curve: points (xyz) in object coordinates of the mesh, kind
    // (0 band, 1 edge loop, 2 crease, 3 off), on/off, width (< 0 automatic) and
    // strength 0..100.
    void addCurve(std::span<const float> pointsXyz, int kind, bool active, float width, float strength);

    [[nodiscard]] bool hasMesh() const noexcept { return hasMesh_; }
    [[nodiscard]] int pointCount() const noexcept { return static_cast<int>(points_.size() / 3); }
    [[nodiscard]] int faceCount() const noexcept { return static_cast<int>(faceSizes_.size()); }
    [[nodiscard]] int cornerCount() const noexcept { return static_cast<int>(cornerPoints_.size()); }

    [[nodiscard]] const std::vector<float>& points() const noexcept { return points_; }
    [[nodiscard]] const std::vector<std::int32_t>& faceSizes() const noexcept { return faceSizes_; }
    [[nodiscard]] const std::vector<std::int32_t>& cornerPoints() const noexcept { return cornerPoints_; }
    [[nodiscard]] const std::optional<std::vector<std::int32_t>>& faceValues(arm_face_values kind) const;
    [[nodiscard]] const std::optional<std::vector<std::int32_t>>& sharpEdges() const noexcept { return sharpEdges_; }
    [[nodiscard]] const std::optional<std::vector<float>>& cornerNormals() const noexcept { return cornerNormals_; }
    [[nodiscard]] bool cornerNormalsShared() const noexcept { return cornerNormalsShared_; }
    [[nodiscard]] const std::optional<std::vector<float>>& cornerColors() const noexcept { return cornerColors_; }
    [[nodiscard]] const std::optional<std::vector<float>>& pointDensity() const noexcept { return pointDensity_; }

    struct Uv {
        std::string name;
        std::vector<float> values;
    };
    [[nodiscard]] const std::vector<Uv>& uvSets() const noexcept { return uvSets_; }

    struct Curve {
        std::vector<float> points;
        int kind = ARM_CURVE_EDGE_LOOP;
        bool active = true;
        float width = -1.0f;
        float strength = 50.0f;
    };
    [[nodiscard]] const std::vector<Curve>& curves() const noexcept { return curves_; }

private:
    void requireMesh() const;

    bool hasMesh_ = false;
    std::vector<float> points_;
    std::vector<std::int32_t> faceSizes_;
    std::vector<std::int32_t> cornerPoints_;
    std::array<std::optional<std::vector<std::int32_t>>, 3> faceValues_;
    std::optional<std::vector<std::int32_t>> sharpEdges_;
    std::optional<std::vector<float>> cornerNormals_;
    bool cornerNormalsShared_ = true;
    std::optional<std::vector<float>> cornerColors_;
    std::optional<std::vector<float>> pointDensity_;
    std::vector<Uv> uvSets_;
    std::vector<Curve> curves_;
};

} // namespace areko::remesher::input
