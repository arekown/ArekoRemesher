#pragma once

namespace areko::remesher::field {

// Radii of the spheres the curvature is measured in: the smallest and
// largest radius and the number of levels between them.
class CurvatureRadii {
public:
    // From the part's target edge length L and area A (float):
    //  - without border data L is at least sqrtf(A / 10000) (maxss),
    //  - r_min = 0.33 L, r_max = 1.75 L; r_min at most 0.8 r_max (in double),
    //  - with several parts both radii are clamped to 0.5 sqrtf(A / 24); the
    //    ratio q = r_max / r_min of before the clamp decides the levels, and
    //    when both radii were clamped r_min becomes r_max / q,
    //  - q == 1 gives one level, otherwise 2 (q < 1.7), 3 (< 3), 4 (< 4),
    //    5 (< 7) or 6.
    static CurvatureRadii compute(float edgeLength, float area, bool borderData, bool severalParts);

    [[nodiscard]] float minRadius() const noexcept { return rMin_; }
    [[nodiscard]] float maxRadius() const noexcept { return rMax_; }
    [[nodiscard]] int levelCount() const noexcept { return levels_; }

    // Radius of level l: the mean of the linear and the geometric step from
    // r_min to r_max (float, log10f/powf); r_max with a single level.
    [[nodiscard]] float levelRadius(int level) const;
    // Smallest curvature that counts on a sphere of radius r: 1 / (r * 3.7978768).
    [[nodiscard]] static float minCurvature(float radius);

private:
    float rMin_ = 0.0f;
    float rMax_ = 0.0f;
    int levels_ = 1;
};

} // namespace areko::remesher::field
