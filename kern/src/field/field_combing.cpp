#include "field_combing.h"

#include "axis_rotation.h"

#include "util/index_queue.h"
#include "util/rounding.h"

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

namespace {

// The first frame axis is turned about the normal with the products and the
// sums in float; axes shorter than 1e-8 give no rotation.
struct FrameTurn {
    using Real = float;
    static constexpr double shortestAxis = 1e-8;
    static constexpr bool addsZero = false;
};

// k modulo 4 in 0 .. 3.
int quarterOf(int k)
{
    return ((k % 4) + 4) % 4;
}

} // namespace

void FieldCombing::comb(const Mesh& mesh, const CrossField& field, const std::vector<int>* faceMirror, CutGraph& out)
{
    const Adjacency& adjacency = mesh.adjacency();
    const std::size_t faceCount = static_cast<std::size_t>(mesh.faceCount());
    out.rotation.assign(faceCount, 0);
    BitArray visited(faceCount);
    util::IndexQueue queue(faceCount);
    for (std::size_t start = 0; start < faceCount; ++start) {
        if (visited.test(start)) {
            continue;
        }
        visited.set(start);
        queue.pushBack(static_cast<int>(start));
        while (!queue.empty()) {
            const int t = queue.pop();
            for (int c = 3 * t; c < 3 * t + 3; ++c) {
                const int e = adjacency.cornerEdge(c);
                if (e < 0 || out.cutEdgeBits.testBounded(static_cast<std::size_t>(e))) {
                    continue;
                }
                const auto faces = adjacency.innerEdgeFaces(e);
                if (!faces) {
                    continue;
                }
                const bool first = (*faces)[0] == t;
                const std::size_t n = static_cast<std::size_t>(first ? (*faces)[1] : (*faces)[0]);
                if (visited.test(n)) {
                    continue;
                }
                const int jump = field.periodJump[static_cast<std::size_t>(e)];
                const float turns = static_cast<float>(out.rotation[static_cast<std::size_t>(t)]) +
                                    static_cast<float>(first ? -jump : jump);
                out.rotation[n] = quarterOf(util::roundHalfAway(turns));
                visited.set(n);
                queue.pushBack(static_cast<int>(n));
            }
        }
    }
    if (faceMirror != nullptr) {
        out.mirrorHalf.forEachSet([&](std::size_t t) {
            // Every triangle of the half has a mirror triangle.
            const int mirrored = (*faceMirror)[t];
            if (mirrored >= 0) {
                out.rotation[t] = 2 - out.rotation[static_cast<std::size_t>(mirrored)];
            }
        });
    }

    const std::vector<Vec3f>& normals = mesh.triangleGeometry().normals();
    out.combedU.resize(faceCount);
    out.combedV.resize(faceCount);
    for (std::size_t t = 0; t < faceCount; ++t) {
        const float angle = static_cast<float>(out.rotation[t]) * quarterTurnF + field.angle[t];
        const Vec3f e1 = field.frameE1[t];
        const Vec3f turned = AxisRotation<FrameTurn>(normals[t], angle).apply(e1);
        const float u = dot(e1, turned);
        const float v = dot(turned, field.frameE2[t]);
        out.combedU[t] = {u, v};
        out.combedV[t] = {-v, u};
    }
}

void FieldCombing::residualJumps(const Adjacency& adjacency, const CrossField& field, CutGraph& out)
{
    const std::size_t edgeCount = static_cast<std::size_t>(adjacency.edgeCount());
    out.residualJump.assign(edgeCount, 0);
    for (std::size_t e = 0; e < edgeCount; ++e) {
        if (const auto faces = adjacency.innerEdgeFaces(static_cast<int>(e))) {
            out.residualJump[e] = out.rotation[static_cast<std::size_t>((*faces)[1])] -
                                  out.rotation[static_cast<std::size_t>((*faces)[0])] + field.periodJump[e];
        }
    }
}

} // namespace areko::remesher::field
