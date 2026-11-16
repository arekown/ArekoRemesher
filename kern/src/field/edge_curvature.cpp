#include "edge_curvature.h"

#include "mesh/measures.h"

#include <cmath>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

namespace {

// 20 degrees in double, the limit of a sharp crease edge.
constexpr double sharpAngle = 0.3490658503988659;

// Unsigned dihedral angle from the cosine between the normals (angleOfCosine,
// so a NaN cosine stays NaN); 0 for a zero denominator.
float unsignedAngle(double cosine, bool valid)
{
    return valid ? static_cast<float>(angleOfCosine(cosine)) : 0.0f;
}

} // namespace

EdgeCurvature::EdgeCurvature(const Mesh& mesh)
    : lengths_(edgeLengths(mesh))
{
    const Adjacency& adjacency = mesh.adjacency();
    const std::vector<Vec3f>& points = mesh.points();
    const std::vector<Vec3f>& normals = mesh.triangleGeometry().normals();
    const int edgeCount = adjacency.edgeCount();
    angles_.assign(static_cast<std::size_t>(edgeCount), 0.0f);
    directions_.assign(static_cast<std::size_t>(edgeCount), Vec3f{});
    tensors_.assign(6 * static_cast<std::size_t>(edgeCount), 0.0);
    for (int e = 0; e < edgeCount; ++e) {
        const std::size_t ei = static_cast<std::size_t>(e);
        const auto ends = adjacency.edgeVertices(e);
        const float l = lengths_[ei];
        Vec3f dir;
        if (l != 0.0f) {
            const Vec3f a = points[static_cast<std::size_t>(ends[0])];
            const Vec3f b = points[static_cast<std::size_t>(ends[1])];
            dir = {(b.x - a.x) / l, (b.y - a.y) / l, (b.z - a.z) / l};
        }
        directions_[ei] = dir;
        float angle = 0.0f;
        const auto corners = adjacency.edgeCorners(e);
        if (corners.size() == 2) {
            const int c0 = corners[0];
            const int c1 = corners[1];
            const int f0 = c0 / 3;
            const int f1 = c1 / 3;
            const Vec3d n0 = wide(normals[static_cast<std::size_t>(f0)]);
            const Vec3d n1 = wide(normals[static_cast<std::size_t>(f1)]);
            const Vec3f o0 = points[static_cast<std::size_t>(adjacency.vertex(3 * f0 + (c0 + 2) % 3))];
            const Vec3f o1 = points[static_cast<std::size_t>(adjacency.vertex(3 * f1 + (c1 + 2) % 3))];
            const double side = dot(wide(o1) - wide(o0), n0);
            const double denominator = length(n0) * length(n1);
            const bool valid = denominator != 0.0;
            const double cosine = valid ? dot(n0, n1) / denominator : 0.0;
            angle = side < 0.0 ? unsignedAngle(cosine, valid) : -unsignedAngle(cosine, valid);
        }
        angles_[ei] = angle;
        const double w = static_cast<double>(l * angle);
        const double x = static_cast<double>(dir.x);
        const double y = static_cast<double>(dir.y);
        const double z = static_cast<double>(dir.z);
        const double xw = x * w;
        double* t = tensors_.data() + 6 * ei;
        t[0] = x * xw;
        t[1] = xw * y;
        t[2] = xw * z;
        t[3] = (y * w) * y;
        t[4] = (w * y) * z;
        t[5] = (w * z) * z;
    }
}

BitArray EdgeCurvature::sharpPoints(const Mesh& mesh, const std::vector<curves::GuideCurve>& curves) const
{
    const Adjacency& adjacency = mesh.adjacency();
    BitArray sharp(static_cast<std::size_t>(mesh.pointCount()));
    for (const curves::GuideCurve& curve : curves) {
        if (curve.kind != ARM_CURVE_CREASE || !curve.active) {
            continue;
        }
        const std::vector<int>& chain = curve.cutVertices;
        for (std::size_t i = 0; i + 1 < chain.size(); ++i) {
            const int v = chain[i];
            const int w = chain[i + 1];
            const int e = adjacency.edgeBetween(v, w);
            if (e != -1 && !(sharpAngle > wide(std::fabs(angles_[static_cast<std::size_t>(e)])))) {
                sharp.setBounded(static_cast<std::size_t>(v));
                sharp.setBounded(static_cast<std::size_t>(w));
            }
        }
    }
    return sharp;
}

} // namespace areko::remesher::field
