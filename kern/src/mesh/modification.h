#pragma once

#include "attribute_maps.h"
#include "bit_array.h"
#include "csr.h"
#include "index_map.h"
#include "topology.h"
#include "vec3.h"

#include <optional>
#include <vector>

namespace areko::remesher::mesh {

class Mesh;

// Description of one mesh change. Every stage changes meshes only through it,
// so maps, UV sets, normals and sharp edges follow in one place.
struct Modification {
    // Points: complete replacement (receives the old points) or points to append.
    std::optional<std::vector<Vec3f>> replacePoints;
    std::vector<Vec3f> appendPoints;
    // Interpolation sources (old point numbers) per appended point, optional weights.
    std::optional<Csr> pointSources;
    std::optional<std::vector<float>> pointWeights;

    // Faces: complete corner list and/or face starts (receive the old lists),
    // mode afterwards (empty = unchanged).
    std::optional<std::vector<int>> cornerVertex;
    std::optional<std::vector<int>> faceStart;
    std::optional<FaceMode> faceMode;

    // Faces to delete and faces to append. The appended faces refer to the
    // points of the changed mesh and must not carry points of their own.
    std::optional<BitArray> deleteFaces;
    const Mesh* appendFaces = nullptr;

    // Renumbering per element kind; faces.newToOld may name the origin of the
    // appended faces, halfEdges.newToOld the origin of the corners.
    IndexMap points;
    IndexMap faces;
    IndexMap corners;
    IndexMap halfEdges;
    IndexMap edges;
    IndexMap uvEdges;

    // Further maps of the caller that follow the same change.
    AttributeMaps* extraMaps = nullptr;
};

class MeshModifier {
public:
    // Applies the change. Returns false and leaves the mesh untouched when there
    // is nothing to do or the appended faces bring their own points.
    static bool apply(Mesh& mesh, Modification& change);
};

} // namespace areko::remesher::mesh
