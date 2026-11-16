#pragma once

#include "mesh/vec3.h"
#include "util/rounding.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace areko::remesher::field {

// A quarter turn (pi / 2) as float: the step between the four arms of a cross.
inline constexpr float quarterTurnF = 1.5707964f;

// Nearest whole number of quarter turns between the crosses of the two
// triangles of an edge: the angle in the second triangle minus the angle in
// the first plus the connection of the edge, over pi / 2 in double, rounded
// to float and then half away from zero.
[[nodiscard]] inline int nearestQuarterTurns(double first, double second, float connection)
{
    const float turns = static_cast<float>(((second - first) + mesh::wide(connection)) / mesh::halfPiD);
    return util::roundHalfAway(mesh::wide(turns));
}

// Points where the crosses do not close up around the point. `index` is the
// sum of the quarter turns across the edges around the point, modulo 4; the
// crosses themselves turn the other way, so 3 is a quarter turn forward
// (three quads will meet there) and 1 a quarter turn back (five quads).
// `fullTurn` marks a point around which the crosses make one whole turn,
// which the count modulo 4 cannot see.
struct Singularities {
    std::vector<int> points;            // ascending
    std::vector<int> index;
    std::vector<std::uint8_t> fullTurn;
    int indexSum = 0;

    // Turn of the crosses around entry i in quarter turns: 1 three quads,
    // -1 five; an index of 2 (half a turn either way) counts as 2, a whole
    // turn adds 4. Over a closed surface these add up to four times its Euler
    // characteristic.
    [[nodiscard]] int quarterTurns(std::size_t i) const
    {
        const int turn = index[i] == 3 ? 1 : (index[i] == 1 ? -1 : index[i]);
        return turn + (fullTurn[i] != 0 ? 4 : 0);
    }
};

// Cross field of one part: per triangle a frame in its plane and the angle of
// the cross against the first frame axis, per edge the number of quarter
// turns between the crosses of its two triangles.
struct CrossField {
    std::vector<mesh::Vec3f> frameE1;       // along the reference edge of the triangle
    std::vector<mesh::Vec3f> frameE2;       // normal x frameE1
    std::vector<std::uint8_t> startCorner;  // corner the reference edge starts at (255: none of the three)
    std::vector<float> referenceLength;     // length of the reference edge
    std::vector<double> targetAngle;        // angle of the curvature direction in the frame
    std::vector<float> angle;               // angle of the cross in the frame
    std::vector<int> periodJump;            // per edge
    Singularities singularities;

    // The vector with the given two coordinates in the frame of triangle t:
    // per axis the part along the second frame axis plus the part along the first.
    [[nodiscard]] mesh::Vec3f inSpace(std::size_t t, const std::array<float, 2>& coordinates) const
    {
        return frameE2[t] * coordinates[1] + frameE1[t] * coordinates[0];
    }

    // One arm of the cross of triangle t (unit length); the others follow by
    // quarter turns about the normal.
    [[nodiscard]] mesh::Vec3f direction(std::size_t t) const
    {
        const float c = std::cos(angle[t]);
        const float s = std::sin(angle[t]);
        return inSpace(t, {c, s});
    }
};

} // namespace areko::remesher::field
