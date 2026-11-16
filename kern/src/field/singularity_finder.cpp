#include "singularity_finder.h"

#include "mesh/measures.h"

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

namespace {

// Field and surface turn alike when they differ by less than this.
constexpr float sameTurn = 0.05f;

// Remainder modulo 4 in 0 .. 3, also for negative numbers.
int modulo4(int value) { return static_cast<int>(static_cast<unsigned>(value) & 3u); }

} // namespace

Singularities SingularityFinder::find(const Mesh& mesh, const CrossField& field, const std::vector<float>& connection,
                                      const prep::SymmetryCache* mirror, bool mirrorValid)
{
    const Adjacency& adjacency = mesh.adjacency();
    const std::vector<Vec3f>& points = mesh.points();
    const std::vector<Vec3f>& normals = mesh.triangleGeometry().normals();
    const ConnectedParts& parts = mesh.connectedParts();
    Singularities out;
    for (int v = 0; v < mesh.pointCount(); ++v) {
        const auto fan = adjacency.pointCorners(v);
        if (fan.empty()) {
            // A point without triangles has nothing to turn around.
            continue;
        }
        double surfaceTurn = 0.0;
        double fieldTurn = 0.0;
        int turns = 0;
        bool border = false;
        for (const int c : fan) {
            const int t = c / 3;
            const int before = adjacency.prev(c);
            const Vec3f p = points[static_cast<std::size_t>(adjacency.vertex(c))];
            const Vec3f a = points[static_cast<std::size_t>(adjacency.vertex(before))] - p;
            const Vec3f b = points[static_cast<std::size_t>(adjacency.vertex(adjacency.next(c)))] - p;
            const Vec3f n = normals[static_cast<std::size_t>(t)];
            if (!(n == Vec3f{})) {
                // The corner angle from b to a about the normal, in double.
                surfaceTurn = surfaceTurn + signedVectorAngleAsDouble(b, a, cross(b, a), n);
            }
            // The edge from the previous corner to this one ends at v.
            const int e = adjacency.cornerEdge(before);
            const auto faces = e >= 0 ? adjacency.innerEdgeFaces(e) : std::nullopt;
            if (!faces) {
                border = true;
                continue;
            }
            const auto [f0, f1] = *faces;
            const int jump = field.periodJump[static_cast<std::size_t>(e)];
            int quarter = modulo4(jump);
            if (quarter == 3) {
                quarter = -1;
            }
            const double across =
                wide(((connection[static_cast<std::size_t>(e)] - field.angle[static_cast<std::size_t>(f0)]) +
                      field.angle[static_cast<std::size_t>(f1)]) -
                     static_cast<float>(jump) * quarterTurnF);
            if (t == f0) {
                turns += quarter;
                fieldTurn = fieldTurn + across;
            } else {
                turns -= quarter;
                fieldTurn = fieldTurn - across;
            }
        }
        if (border) {
            continue;
        }
        const int index = modulo4(turns);
        const float difference = static_cast<float>(fieldTurn) - static_cast<float>(surfaceTurn);
        bool fullTurn = difference > -sameTurn && sameTurn > difference;
        if (fullTurn && mirror != nullptr && parts.partOfPoint(v) != -1 && mirrorValid) {
            const std::size_t i = static_cast<std::size_t>(v);
            const std::vector<int>& pairs = mirror->pointMirror();
            if (!mirror->pointSide().testBounded(i) && (i < pairs.size() ? pairs[i] : -1) != v) {
                fullTurn = false;
            }
        }
        if (fullTurn || index != 0) {
            out.points.push_back(v);
            out.index.push_back(index);
            out.fullTurn.push_back(fullTurn ? 1 : 0);
            out.indexSum += index;
        }
    }
    return out;
}

} // namespace areko::remesher::field
