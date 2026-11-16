#include "cross_frames.h"

#include "axis_rotation.h"

#include "mesh/measures.h"

#include <cmath>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

namespace {

using Rotation = AxisRotation<PlaneTurn>;

// Turns the plane with normal `from` into the plane with normal `to`.
Rotation planeRotation(Vec3f from, Vec3f to)
{
    const Vec3f crossed = cross(from, to);
    if (!(wide(dot(crossed, crossed)) < 1e-7)) {
        Vec3f axis = crossed;
        normalizeInPlace(axis);
        return Rotation(axis, signedVectorAngle(from, to, crossed, axis));
    }
    if (dot(to, from) < 0.0f) {
        // Opposite normals: half a turn about any axis in the plane.
        return Rotation(cross(smallestAxis(from), from), piF);
    }
    return Rotation();
}

} // namespace

void CrossFrames::compute(const Mesh& mesh, prep::SymmetryCache* mirror, CrossField& field)
{
    const Topology& topology = mesh.topology();
    const Adjacency& adjacency = mesh.adjacency();
    const std::vector<Vec3f>& points = mesh.points();
    const std::vector<Vec3f>& normals = mesh.triangleGeometry().normals();
    const std::size_t faceCount = static_cast<std::size_t>(mesh.faceCount());
    field.frameE1.assign(faceCount, Vec3f{});
    field.frameE2.assign(faceCount, Vec3f{});
    field.startCorner.assign(faceCount, 0);
    field.referenceLength.assign(faceCount, 0.0f);
    const std::vector<int>* faceMirror = mirror != nullptr ? &mirror->faceMirror(mesh) : nullptr;
    const std::vector<int>* edgeMirror = mirror != nullptr ? &mirror->edgeMirror(adjacency) : nullptr;
    auto mirrorOfEdge = [&](int e) {
        return e >= 0 && static_cast<std::size_t>(e) < edgeMirror->size() ? (*edgeMirror)[static_cast<std::size_t>(e)] : -2;
    };
    for (std::size_t t = 0; t < faceCount; ++t) {
        const int first = 3 * static_cast<int>(t);
        int from = 0;
        int to = 1;
        if (mirror != nullptr && t < faceMirror->size() && (*faceMirror)[t] != -1 &&
            mirror->sideOfFace(topology, static_cast<int>(t)) == 0) {
            // Mirror image of the partner's reference edge.
            const int image = mirrorOfEdge(adjacency.cornerEdge(3 * (*faceMirror)[t]));
            if (image == adjacency.cornerEdge(first)) {
                field.startCorner[t] = 0;
            } else if (image == adjacency.cornerEdge(first + 1)) {
                field.startCorner[t] = 1;
                from = 1;
                to = 2;
            } else {
                // Also taken when the image is no edge of this triangle.
                field.startCorner[t] = image == adjacency.cornerEdge(first + 2) ? 2 : 255;
                from = 2;
                to = 0;
            }
        }
        Vec3f e1 = points[static_cast<std::size_t>(topology.vertex(first + to))] -
                   points[static_cast<std::size_t>(topology.vertex(first + from))];
        field.referenceLength[t] = normalizeWithLength(e1);
        Vec3f e2 = cross(normals[t], e1);
        normalizeInPlace(e2);
        field.frameE1[t] = e1;
        field.frameE2[t] = e2;
    }
}

std::vector<float> CrossFrames::connection(const Mesh& mesh, const std::vector<Vec3f>& frameE1)
{
    const Adjacency& adjacency = mesh.adjacency();
    const std::vector<Vec3f>& normals = mesh.triangleGeometry().normals();
    std::vector<float> out(static_cast<std::size_t>(adjacency.edgeCount()), 0.0f);
    for (int e = 0; e < adjacency.edgeCount(); ++e) {
        const auto faces = adjacency.innerEdgeFaces(e);
        if (!faces) {
            continue;
        }
        const std::size_t f = static_cast<std::size_t>((*faces)[0]);
        const std::size_t g = static_cast<std::size_t>((*faces)[1]);
        const Vec3f moved = planeRotation(normals[g], normals[f]).apply(frameE1[g]);
        float angle = signedVectorAngle(frameE1[f], moved, cross(frameE1[f], moved), normals[f]);
        if (angle > piF) {
            angle = angle - fullTurnF;
        }
        out[static_cast<std::size_t>(e)] = angle;
    }
    return out;
}

} // namespace areko::remesher::field
