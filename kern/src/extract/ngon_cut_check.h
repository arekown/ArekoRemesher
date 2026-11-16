#pragma once

#include "ngon_new_points.h"
#include "ngon_polygon_list.h"

namespace areko::remesher::extract {

// Closest approach of a ray (origin o, unit direction d) to the segment a-b:
// distance, ray parameter t and segment parameter u in [0, 1].
struct Approach {
    float distance = 0.0f;
    float along = 0.0f;
    float segment = 0.0f;
};

// Whether a chord between two corners of a polygon may split it, and what it costs.
class CutCheck {
public:
    [[nodiscard]] static Approach closest(mesh::Vec3f o, mesh::Vec3f d, mesh::Vec3f a, mesh::Vec3f b);

    struct Chord {
        bool blocked = false;
        float length = 0.0f;
    };
    // Chord from position i to position j of polygon p: blocked when an edge
    // not touching i or j passes close to it inside its span.
    [[nodiscard]] static Chord chord(const NgonPolygonList& list, const NgonNewPoints& points, int p, int i, int j);
    // Multiplies `cost` by (1 + 3 |w1 - w2| / 2 / (pi/4)) at i and at j, where
    // w1, w2 are the angles between the chord and the two edges at that end.
    [[nodiscard]] static float angleFactor(const NgonPolygonList& list, const NgonNewPoints& points, int p, int i, int j,
                                           float length, float cost);
};

} // namespace areko::remesher::extract
