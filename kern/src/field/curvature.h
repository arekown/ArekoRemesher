#pragma once

#include "mesh/bit_array.h"
#include "mesh/vec3.h"

#include <cstdint>
#include <vector>

namespace areko::remesher::field {

// Marks of a triangle in the curvature result: bit 0 the direction is set,
// bit 1 the curvature values are reliable, bits 2/3 the kind of a fixed
// direction (one value, assigned, never combined).
namespace curvature_flags {
inline constexpr std::uint8_t direction = 1;
inline constexpr std::uint8_t reliable = 2;
inline constexpr std::uint8_t fixedMask = 0xc;
inline constexpr std::uint8_t fixedCurve = 5;    // guide curve (with the direction bit)
inline constexpr std::uint8_t fixedBorder = 9;   // border edge
inline constexpr std::uint8_t fixedPlane = 13;   // edge on the mirror plane
} // namespace curvature_flags

// Distance of the points to the guide curves: per point its distance and the
// seed it came from (-1 none), per seed its curve (index in the part's curve
// list) and the point on that curve; per triangle the mean of its three
// distances, each over the width of its curve (FLT_MAX unless all three
// points have a seed).
struct CurveDistance {
    std::vector<float> pointDistance;
    std::vector<int> pointSeed;
    std::vector<int> seedCurve;
    std::vector<int> seedPoint;
    std::vector<float> faceDistance;
};

// Curvature of one part, per triangle: the direction of least curvature
// (zero where none was found; the half that is not computed keeps only its
// fixed directions), the two curvature values (the second with its sign),
// the marks and the trust weight of the direction.
struct Curvature {
    std::vector<mesh::Vec3f> direction;
    std::vector<float> lambda1;
    std::vector<float> lambda2;
    std::vector<std::uint8_t> flags;
    std::vector<float> weight;

    // Triangles of the mirror half that are not computed and how many they are;
    // points whose tensor is not computed (mirror points of that half).
    mesh::BitArray mirrorFaces;
    int mirrorFaceCount = 0;
    mesh::BitArray mirrorPoints;
    // Points that end the flooding of the vertex tensors (crease curves).
    mesh::BitArray sharpPoints;
    CurveDistance curveDistance;
};

} // namespace areko::remesher::field
