#pragma once

#include "triangulator.h"

#include "mesh/bit_array.h"
#include "mesh/mesh.h"

#include <optional>
#include <utility>
#include <vector>

namespace areko::remesher::prep {

// Mirror pairs of a mirror-symmetric mesh, from a cut at a plane or from the
// detection:
//  - per point its mirror point (itself on the plane, -1 none) and its side,
//  - per face its mirror face and its side,
//  - per edge its mirror edge and the edges on the plane, built from the
//    point pairs on first use when the cache does not bring them.
// The pairs refer to the mesh they were made for; a later change of points or
// faces makes them invalid, and a new detection takes their place.
class SymmetryCache {
public:
    struct Maps {
        std::vector<int> points;
        mesh::BitArray pointSide; // cut: 0 computed half, 1 mirrored half; detection: 1 towards the normal
        std::vector<int> faces;
        mesh::BitArray faceSide;
    };

    // Plane fixed by the cut: a new detection searches only along its axis
    // and keeps normal and point.
    struct LockedPlane {
        int axis = -1;
        mesh::Vec3f normal;
        mesh::Vec3f point;
    };

    // Pairs of a cut: the face pairs count as not yet checked, the edge pairs
    // follow from the point pairs.
    SymmetryCache(Maps maps, const LockedPlane& plane) : maps_(std::move(maps)), plane_(plane) {}

    // Pairs of a detection with its edge pairs and the number of paired faces.
    SymmetryCache(Maps maps, std::vector<int> edgeMirror, int mirroredFaces, std::optional<LockedPlane> plane)
        : maps_(std::move(maps)), edgeMirror_(std::move(edgeMirror)), edgeMirrorBuilt_(true),
          mirroredFaces_(mirroredFaces), plane_(plane)
    {
    }

    // Point pairs and sides only (face pairs not yet checked, edge pairs from
    // the points), with the plane given; for copies onto other meshes.
    [[nodiscard]] static SymmetryCache fromPoints(std::vector<int> points, mesh::BitArray pointSide,
                                                  const std::optional<LockedPlane>& plane);

    [[nodiscard]] const std::vector<int>& pointMirror() const noexcept { return maps_.points; }
    [[nodiscard]] const mesh::BitArray& pointSide() const noexcept { return maps_.pointSide; }
    // The face pairs as they are stored.
    [[nodiscard]] const std::vector<int>& faceMirror() const noexcept { return maps_.faces; }
    [[nodiscard]] const mesh::BitArray& faceSide() const noexcept { return maps_.faceSide; }
    [[nodiscard]] const std::optional<LockedPlane>& lockedPlane() const noexcept { return plane_; }
    // Number of paired faces (a face on itself counts once); -1 not counted.
    [[nodiscard]] int mirroredFaceCount() const noexcept { return mirroredFaces_; }

    // Point and face pairs for the mirror path of the triangulation.
    [[nodiscard]] MirrorLists mirrorLists() const noexcept { return {maps_.points, maps_.faces}; }

    // Edge pairs; without them edge (a, b) pairs with the first edge
    // (m(a), m(b)), edges with a pair keep it.
    [[nodiscard]] const std::vector<int>& edgeMirror(const mesh::Adjacency& adjacency);
    // Edges that mirror onto themselves, ascending.
    [[nodiscard]] const std::vector<int>& planeEdges(const mesh::Adjacency& adjacency);
    // Face pairs; uncounted ones are rebuilt from the point pairs (triangle
    // meshes: the first triangle at edge (m(a), m(b)) with corner m(c), each
    // target face once; other meshes get none).
    [[nodiscard]] const std::vector<int>& faceMirror(const mesh::Mesh& mesh);
    // Whether every face has a mirror face (rebuilding the face pairs if needed).
    [[nodiscard]] bool allFacesMirrored(const mesh::Mesh& mesh);

    // Side of face f from its corners off the plane: 1 all on side 1, 0 all
    // on side 0, 2 as many on both, -1 otherwise or without point pairs.
    [[nodiscard]] int sideOfFace(const mesh::Topology& topology, int f) const;

private:
    void buildEdgeMirror(const mesh::Adjacency& adjacency);

    Maps maps_;
    std::vector<int> edgeMirror_;
    bool edgeMirrorBuilt_ = false;
    int mirroredFaces_ = -1;
    std::optional<std::vector<int>> planeEdges_;
    std::optional<LockedPlane> plane_;
};

} // namespace areko::remesher::prep
