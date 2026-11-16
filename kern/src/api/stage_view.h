#pragma once

#include "areko_remesher.h"
#include "mesh/mesh.h"
#include "prep/net_lines.h"

#include <array>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace areko::remesher::api {

// Flat copies of the data of one computed stage, in the shape the C
// interface hands out.
class StageView {
public:
    // Lines that stand for sharp edges are shown as sharp edges; from the
    // stage "feature lines" on, every line is handed out as well. `mirrorSides`
    // (one value per face or empty) becomes the face values of the mirror side.
    StageView(arm_stage stage, const mesh::Mesh& mesh, std::span<const prep::NetLine> lines,
              std::span<const int> mirrorSides);

    // Part number per face, budget and skip mark per part (pairs), kind per line.
    void setParts(std::span<const int> faceParts, std::span<const int> partValues);
    void setLineKinds(std::span<const int> kinds);
    // Curvature direction per face (xyz).
    void setFaceDirections(std::span<const float> directions);
    // Cross field: one arm of the cross per face (xyz), and the singular
    // points as pairs (point, turn of the crosses in quarter turns).
    void setCrossField(std::span<const float> directions, std::span<const int> singularities);
    // Cut edges as point pairs.
    void setCutEdges(std::span<const int> pointPairs);
    // Size field: two target sizes per face, four values per part.
    void setSizeField(std::span<const float> faceSizes, std::span<const float> partSizes);
    // Grid coordinates: two per corner, and per part 1 when it has any.
    void setGrid(std::span<const float> cornerCoordinates, std::span<const int> partSolved);

    [[nodiscard]] arm_stage stage() const noexcept { return stage_; }
    [[nodiscard]] arm_mesh_view meshView() const noexcept;
    [[nodiscard]] const std::optional<std::vector<std::int32_t>>& faceValues(arm_face_values kind) const;
    [[nodiscard]] const std::vector<std::int32_t>& sharpEdgePairs() const noexcept { return sharpEdgePairs_; }
    // Line i runs over the points linePoints()[lineStarts()[i] .. lineStarts()[i + 1]).
    [[nodiscard]] const std::vector<std::int32_t>& lineStarts() const noexcept { return lineStarts_; }
    [[nodiscard]] const std::vector<std::int32_t>& linePoints() const noexcept { return linePoints_; }
    [[nodiscard]] const std::vector<std::int32_t>& lineKinds() const noexcept { return lineKinds_; }
    [[nodiscard]] const std::vector<std::int32_t>& partValues() const noexcept { return partValues_; }
    [[nodiscard]] const std::vector<std::int32_t>& cornerNormal() const noexcept { return cornerNormal_; }
    [[nodiscard]] const std::vector<float>& normals() const noexcept { return normals_; }
    [[nodiscard]] const std::vector<std::uint8_t>& colors() const noexcept { return colors_; }
    [[nodiscard]] const std::vector<float>& density() const noexcept { return density_; }
    [[nodiscard]] const std::vector<float>& faceDirections() const noexcept { return faceDirections_; }
    [[nodiscard]] const std::vector<float>& crossDirections() const noexcept { return crossDirections_; }
    [[nodiscard]] const std::vector<std::int32_t>& singularities() const noexcept { return singularities_; }
    [[nodiscard]] const std::vector<std::int32_t>& cutEdgePairs() const noexcept { return cutEdgePairs_; }
    [[nodiscard]] const std::vector<float>& faceSizes() const noexcept { return faceSizes_; }
    [[nodiscard]] const std::vector<float>& partSizes() const noexcept { return partSizes_; }
    [[nodiscard]] const std::vector<float>& gridCoordinates() const noexcept { return gridCoordinates_; }
    [[nodiscard]] const std::vector<std::int32_t>& gridParts() const noexcept { return gridParts_; }

    struct Uv {
        std::string name;
        std::vector<std::int32_t> cornerUv;
        std::vector<float> coords;
    };
    [[nodiscard]] const std::vector<Uv>& uvSets() const noexcept { return uvSets_; }

private:
    arm_stage stage_;
    std::vector<float> points_;
    std::vector<std::int32_t> faceStarts_;
    std::vector<std::int32_t> corners_;
    std::array<std::optional<std::vector<std::int32_t>>, 5> faceValues_;
    std::vector<std::int32_t> sharpEdgePairs_;
    std::vector<std::int32_t> lineStarts_{0};
    std::vector<std::int32_t> linePoints_;
    std::vector<std::int32_t> lineKinds_;
    std::vector<std::int32_t> partValues_;
    std::vector<std::int32_t> cornerNormal_;
    std::vector<float> normals_;
    std::vector<std::uint8_t> colors_;
    std::vector<float> density_;
    std::vector<float> faceDirections_;
    std::vector<float> crossDirections_;
    std::vector<std::int32_t> singularities_;
    std::vector<std::int32_t> cutEdgePairs_;
    std::vector<float> faceSizes_;
    std::vector<float> partSizes_;
    std::vector<float> gridCoordinates_;
    std::vector<std::int32_t> gridParts_;
    std::vector<Uv> uvSets_;
};

} // namespace areko::remesher::api
