#include "settings.h"

#include "api/status.h"
#include "mesh/vec3.h"

#include <algorithm>
#include <cstring>

namespace areko::remesher::settings {

Settings::Settings() noexcept
{
    fillDefaults(block_);
}

void Settings::fillDefaults(arm_settings& block) noexcept
{
    block = arm_settings{};
    block.struct_size = sizeof(arm_settings);
    block.quad_count = 5000;
    block.quad_percentage = 0.0f;
    block.edge_length = 0.0f;
    block.exact_count = 0;
    block.color_density = 0;
    block.attribute_density = 0;
    block.curvature = 50.0f;
    block.points_per_circle = 0.0f;
    block.density_ratio = 0.0f;
    block.max_quad_ratio = 4.0f;
    block.detect_creases = 1;
    block.denoise = 1;
    block.crease_angle = 30.0f;
    block.follow_borders = 1;
    block.border_density_factor = 1.0f;
    block.mirror_at_center = 1;
    block.thread_count = 0;
}

void Settings::assign(const arm_settings& block)
{
    if (block.struct_size < sizeof(std::uint32_t)) {
        throw api::StatusError(ARM_ERROR_INVALID_ARGUMENT, "settings: struct_size too small");
    }
    arm_settings merged;
    fillDefaults(merged);
    const std::size_t n = std::min<std::size_t>(block.struct_size, sizeof(arm_settings));
    std::memcpy(&merged, &block, n);
    merged.struct_size = sizeof(arm_settings);
    block_ = merged;
}

void SettingsConversion::curvatureToDensity(float c, float& pointsPerCircle, float& densityRatio) noexcept
{
    if (c < 0.0f) {
        pointsPerCircle = 4.0f;
        densityRatio = 1.0f;
    } else if (c < 50.0f) {
        pointsPerCircle = (c * 12.0f) / 50.0f + 4.0f;
        densityRatio = (c * 3.0f) / 50.0f + 1.0f;
    } else if (c > 100.0f) {
        pointsPerCircle = 32.0f;
        densityRatio = 12.0f;
    } else {
        pointsPerCircle = ((c - 50.0f) * 16.0f) / 50.0f + 16.0f;
        densityRatio = ((c - 50.0f) * 8.0f) / 50.0f + 4.0f;
    }
}

std::int32_t SettingsConversion::quadCountFromPercentage(int faceCount, float percentage) noexcept
{
    if (!(percentage > 0.0f) || faceCount == 0) {
        return 0;
    }
    return mesh::truncToInt(static_cast<float>(faceCount) * (percentage / 100.0f));
}

float SettingsConversion::degreesToRadians(float degrees) noexcept
{
    return static_cast<float>(static_cast<double>(degrees) * 0.017453292519943295);
}

float SettingsConversion::inverseBorderFactor(float factor) noexcept
{
    return factor > 0.0f ? 1.0f / factor : 1.0f;
}

SymmetrySwitches SettingsConversion::symmetry(const arm_settings& block) noexcept
{
    // Same rules as a text of axis letters: "detect" alone, otherwise the axes,
    // with a leading "keep half" mark.
    SymmetrySwitches s;
    const bool anyAxis = block.mirror_x || block.mirror_y || block.mirror_z;
    if (block.detect_symmetry) {
        s.on = true;
        s.detect = true;
        return s;
    }
    if (!anyAxis && !block.freeze_half) {
        return s;
    }
    s.axis[0] = block.mirror_x != 0;
    s.axis[1] = block.mirror_y != 0;
    s.axis[2] = block.mirror_z != 0;
    if (block.freeze_half) {
        s.halfMode = true;
    } else {
        s.on = true;
        s.mirrorAtEnd = true;
        s.curvatureSymmetry = true;
    }
    return s;
}

DerivedSettings SettingsConversion::derive(const arm_settings& block, int faceCount)
{
    DerivedSettings d;
    d.quadCount = block.quad_count;
    if (block.quad_percentage > 0.0f && faceCount != 0) {
        d.quadCount = quadCountFromPercentage(faceCount, block.quad_percentage);
    }
    d.edgeLength = block.edge_length;
    if (!(d.edgeLength > 0.0f) && d.quadCount <= 0) {
        throw api::StatusError(ARM_ERROR_INVALID_SETTINGS, "neither a quad count nor an edge length is set");
    }
    curvatureToDensity(block.curvature, d.pointsPerCircle, d.densityRatio);
    if (block.points_per_circle > 0.0f) {
        d.pointsPerCircle = block.points_per_circle;
    }
    if (block.density_ratio > 0.0f) {
        d.densityRatio = block.density_ratio;
    }
    d.borderDensityInverse = inverseBorderFactor(block.border_density_factor);
    d.creaseAngleRadians = degreesToRadians(block.crease_angle);
    d.denoiseStrength = block.denoise ? 0.5f : 0.0f;
    d.symmetry = symmetry(block);
    return d;
}

} // namespace areko::remesher::settings
