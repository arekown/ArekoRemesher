#pragma once

#include "areko_remesher.h"
#include "mesh/vec3.h"

#include <cstdint>
#include <vector>

namespace areko::remesher::curves {

// Whether curves of this kind (arm_curve_kind; lines of the preparation carry
// their type in the same numbers) are cut into the mesh.
[[nodiscard]] inline bool cutsIntoMesh(int kind) noexcept
{
    return kind == ARM_CURVE_EDGE_LOOP || kind == ARM_CURVE_CREASE;
}

// Whether curves of this kind lie on the mesh and steer a band beside them.
[[nodiscard]] inline bool isBand(int kind) noexcept
{
    return kind == ARM_CURVE_BAND;
}

// Whether curves of this kind take part in the field steps (all but "off").
[[nodiscard]] inline bool counts(int kind) noexcept
{
    return kind != ARM_CURVE_OFF;
}

// Where a curve point lies on a mesh: on a point, on an edge or inside a face.
struct Location {
    enum Kind : int { None = 0, Point = 1, Edge = 2, Face = 4 };
    int kind = None;
    int index = -1;

    friend bool operator==(const Location&, const Location&) = default;
};

// A curve laid onto a triangle mesh between two located points: per entry the
// edge under it (-1 inside a face), its parameter along the edge (0 first, 1
// second end point) and its position; per piece entry j -> j+1 the face it
// crosses (-1 unknown); per entry an origin code (empty without codes).
struct CurveSegment {
    std::vector<mesh::Vec3f> points;
    std::vector<int> edges;
    std::vector<float> t;
    std::vector<int> faces;
    std::vector<int> sources;

    [[nodiscard]] int size() const noexcept { return static_cast<int>(edges.size()); }
    [[nodiscard]] bool empty() const noexcept { return edges.empty(); }
    void clear()
    {
        points.clear();
        edges.clear();
        t.clear();
        faces.clear();
        sources.clear();
    }
};

// Points with their locations on a mesh (one location per point).
struct LocatedPolyline {
    std::vector<mesh::Vec3f> points;
    std::vector<Location> locations;
    std::vector<int> faces; // face of each point in the mesh it was located on
};

// One guide curve: what the host or a line of the preparation brought, and
// what projecting and cutting made of it on one part.
struct GuideCurve {
    // Host curves: index in the input list; lines of the preparation: -2, on a
    // part -1.
    int index = -1;
    std::vector<mesh::Vec3f> points;
    int kind = ARM_CURVE_EDGE_LOOP;
    bool active = true;
    float strength = 0.5f; // stored only; no step reads it
    float width = -1.0f;   // < 0: automatic
    std::uint8_t lineFlags = 0;
    // Locations of the points (lines of the preparation: their points) and
    // how many of them are points; -1: the points are located anew.
    std::vector<Location> locations;
    int pointLocationCount = -1;

    // On a part: the located points, the pieces on the part mesh after
    // projecting, and the working polyline: for curves cut into the mesh the
    // cut chain (points of the part mesh and their positions), for bands the
    // source points laid onto the part mesh with the triangle under each
    // (laid once, `bandPrepared`).
    LocatedPolyline source;
    bool sourcePrepared = false;
    std::vector<CurveSegment> segments;
    bool projected = false;
    std::vector<int> cutVertices;
    std::vector<mesh::Vec3f> cutPoints;
    std::vector<int> cutFaces;
    bool bandPrepared = false;

    [[nodiscard]] bool fromHost() const noexcept { return index >= 0; }
};

} // namespace areko::remesher::curves
