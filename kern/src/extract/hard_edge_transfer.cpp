#include "hard_edge_transfer.h"

#include "areko_remesher.h"
#include "field/curve_edges.h"
#include "mesh/sharp_edges.h"

#include <cmath>
#include <utility>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

constexpr float startAngle = 10.0f;
constexpr float widestAngle = 0.7853982f;              // 45 degrees, float
constexpr double sharpWinsBelow = 0.10471975511965978; // 6 degrees, double
constexpr float sharpWeight = 0.6f;

// Whether the edge from a towards b follows a sharp edge of the triangle
// mesh: the best candidate after all edges at a is sharp. An edge whose dot
// product is not a number is skipped; kept, it would get a NaN angle like an
// edge whose cosine is not a number, and a NaN angle never becomes the best
// one, so skipping it changes nothing.
bool followsSharpEdge(const Mesh& triangles, const BitArray& sharp, int a, int b)
{
    const Adjacency& adjacency = triangles.adjacency();
    const auto& points = triangles.points();
    const Vec3f pa = points[static_cast<std::size_t>(a)];
    const Vec3f u = points[static_cast<std::size_t>(b)] - pa;
    int best = -1;
    float bestAngle = startAngle;
    for (const int e : adjacency.pointEdges(a)) {
        const bool inSharp = sharp.testBounded(static_cast<std::size_t>(e));
        const Vec3f w = points[static_cast<std::size_t>(adjacency.otherEnd(e, a))] - pa;
        const float along = dot(u, w);
        if (!(0.0f < along)) {
            continue;
        }
        const double lengths = lengthAsDouble(w) * lengthAsDouble(u);
        float angle = 0.0f;
        if (lengths != 0.0) {
            const double c = std::fabs(static_cast<double>(along) / lengths);
            if (!(c >= 1.0)) { // a NaN cosine gives a NaN angle
                angle = static_cast<float>(std::acos(c));
                if (widestAngle < angle) {
                    continue;
                }
            }
        }
        if (inSharp && static_cast<double>(angle) < sharpWinsBelow) {
            best = e;
            bestAngle = angle;
            continue;
        }
        if (inSharp) {
            angle = angle * sharpWeight;
        }
        if (bestAngle > angle) {
            best = e;
            bestAngle = angle;
        }
    }
    // Decided by the final best candidate only.
    return best != -1 && sharp.testBounded(static_cast<std::size_t>(best));
}

} // namespace

void HardEdgeTransfer::run(const Mesh& triangles, std::span<const curves::GuideCurve> curves, RawMesh& raw)
{
    const Adjacency& triangleAdjacency = triangles.adjacency();
    BitArray sharp(static_cast<std::size_t>(triangleAdjacency.edgeCount()));
    bool any = false;
    for (const curves::GuideCurve& curve : curves) {
        const int flags = curve.lineFlags;
        if (curve.kind != ARM_CURVE_CREASE || !curve.active || !((flags & 4) != 0 || (flags & 3) == 1)) {
            continue;
        }
        field::CurveEdges::forEachChainEdge(triangleAdjacency, curve, [&](int e) {
            sharp.set(static_cast<std::size_t>(e));
            any = true;
        });
    }
    if (!any) {
        return;
    }
    const Adjacency& adjacency = raw.mesh.adjacency();
    const std::vector<PointSource>& sources = raw.points.source;
    BitArray result(static_cast<std::size_t>(adjacency.edgeCount()));
    raw.curveEdges.forEachSet([&](std::size_t q) {
        const auto ends = adjacency.edgeVertices(static_cast<int>(q));
        const PointSource& first = sources[static_cast<std::size_t>(ends[0])];
        const PointSource& second = sources[static_cast<std::size_t>(ends[1])];
        if (first.type == 4 || second.type == 4) {
            return;
        }
        bool hard = false;
        if (first.type == 2) {
            hard = sharp.testBounded(static_cast<std::size_t>(first.id));
        } else if (second.type == 2) {
            hard = sharp.testBounded(static_cast<std::size_t>(second.id));
        } else if (first.type == 1 && second.type == 1) {
            hard = followsSharpEdge(triangles, sharp, first.id, second.id);
        }
        if (hard) {
            result.set(q);
        }
    });
    SharpEdges::setEdgeBits(raw.mesh, std::move(result));
}

} // namespace areko::remesher::extract
