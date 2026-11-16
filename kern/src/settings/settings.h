#pragma once

#include "areko_remesher.h"

#include <cstdint>

namespace areko::remesher::settings {

// The settings block of the interface is the single description of every
// setting; this class owns one block and its defaults.
class Settings {
public:
    Settings() noexcept;

    // Defaults of the interface block.
    static void fillDefaults(arm_settings& block) noexcept;
    // Takes the caller's block; fields beyond its struct_size keep the defaults.
    void assign(const arm_settings& block);

    [[nodiscard]] const arm_settings& block() const noexcept { return block_; }

private:
    arm_settings block_;
};

struct SymmetrySwitches {
    bool on = false;          // symmetric computation
    bool detect = false;      // find the symmetry (no axis given)
    bool mirrorAtEnd = false; // mirror the half at the end
    bool halfMode = false;    // keep the given half
    bool axis[3] = {false, false, false};
    bool curvatureSymmetry = false;
};

// Values the stages read, derived from the block once per run.
struct DerivedSettings {
    std::int32_t quadCount = 0;
    float edgeLength = 0.0f;
    float pointsPerCircle = 0.0f;
    float densityRatio = 0.0f;
    float borderDensityInverse = 1.0f;
    float creaseAngleRadians = 0.0f;
    float denoiseStrength = 0.0f;
    SymmetrySwitches symmetry;
};

class SettingsConversion {
public:
    // Curvature 0..100 -> points per circle and density ratio (float):
    //   C < 0: 4, 1;  0 <= C < 50: C*12/50 + 4, C*3/50 + 1;
    //   50 <= C <= 100 (also NaN): (C-50)*16/50 + 16, (C-50)*8/50 + 4;  C > 100: 32, 12.
    static void curvatureToDensity(float curvature, float& pointsPerCircle, float& densityRatio) noexcept;
    // Percentage of the input faces, after reading the input; 0 when not applicable.
    [[nodiscard]] static std::int32_t quadCountFromPercentage(int faceCount, float percentage) noexcept;
    [[nodiscard]] static float degreesToRadians(float degrees) noexcept;
    [[nodiscard]] static float inverseBorderFactor(float factor) noexcept;
    [[nodiscard]] static SymmetrySwitches symmetry(const arm_settings& block) noexcept;

    // All derived values for an input with `faceCount` faces. Throws a
    // StatusError (invalid settings) when neither a quad count nor an edge
    // length remains.
    [[nodiscard]] static DerivedSettings derive(const arm_settings& block, int faceCount);
};

} // namespace areko::remesher::settings
