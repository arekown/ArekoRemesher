#include "triangulator.h"

#include "polygon_triangulator.h"
#include "mesh/sharp_edges.h"

#include <cmath>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr float twoPiF = 6.2831855f;

// Angle between two vectors: lengths from float square sums, cosine in double.
float angleBetween(Vec3f a, Vec3f b)
{
    const double la = lengthAsDouble(a);
    const double lb = lengthAsDouble(b);
    const double lengths = la * lb;
    if (lengths == 0.0) {
        return 0.0f;
    }
    return static_cast<float>(angleOfCosine(static_cast<double>(dot(a, b)) / lengths));
}

// Corner angles of the diagonal rule: the inner angle per corner, and per quad
// corner the angle sum the corner would get when the diagonal starts there.
struct CornerAngles {
    std::vector<float> inner;
    std::vector<float> split;
    std::vector<char> flatterUnsplit; // per point
};

CornerAngles cornerAngles(const Mesh& mesh)
{
    const Topology& topology = mesh.topology();
    const auto& points = mesh.points();
    CornerAngles out;
    out.inner.assign(static_cast<std::size_t>(topology.cornerCount()), 0.0f);
    out.split = out.inner;
    auto at = [&](int c) { return points[static_cast<std::size_t>(topology.vertex(c))]; };
    for (int f = 0; f < topology.faceCount(); ++f) {
        const int b = topology.faceBegin(f);
        const int n = topology.faceSize(f);
        for (int i = 0; i < n; ++i) {
            const int c = b + i;
            const Vec3f p = at(c);
            const Vec3f next = at(b + (i + 1) % n) - p;
            const Vec3f prev = at(b + (i + n - 1) % n) - p;
            const float inner = angleBetween(next, prev);
            out.inner[static_cast<std::size_t>(c)] = inner;
            if (n == 4) {
                const Vec3f opposite = at(b + (i + 2) % 4) - p;
                out.split[static_cast<std::size_t>(c)] = angleBetween(opposite, prev) + angleBetween(next, opposite);
            } else {
                out.split[static_cast<std::size_t>(c)] = inner;
            }
        }
    }
    const Adjacency& adjacency = mesh.adjacency();
    out.flatterUnsplit.assign(static_cast<std::size_t>(mesh.pointCount()), 0);
    for (int v = 0; v < mesh.pointCount(); ++v) {
        float s1 = 0.0f;
        float s2 = 0.0f;
        for (const int c : adjacency.pointCorners(v)) {
            const float a = out.inner[static_cast<std::size_t>(c)];
            const float s = out.split[static_cast<std::size_t>(c)];
            s1 += a;
            s2 += a < s ? s : a;
        }
        out.flatterUnsplit[static_cast<std::size_t>(v)] = std::fabs(twoPiF - s1) < std::fabs(twoPiF - s2) ? 1 : 0;
    }
    return out;
}

// First corner of the face whose point has exactly two corners and two edges, or -1.
int twoCornerStart(const Adjacency& adjacency, const Topology& topology, int f)
{
    for (int i = 0; i < topology.faceSize(f); ++i) {
        const int v = topology.vertex(topology.faceBegin(f) + i);
        if (adjacency.pointCorners(v).size() == 2 && adjacency.pointEdges(v).size() == 2) {
            return i;
        }
    }
    return -1;
}

int quadStart(const Mesh& mesh, const CornerAngles* angles, int f)
{
    const Topology& topology = mesh.topology();
    const int b = topology.faceBegin(f);
    const int two = twoCornerStart(mesh.adjacency(), topology, f);
    if (two >= 0) {
        return two;
    }
    if (angles != nullptr) {
        auto splitAt = [&](int i) { return angles->split[static_cast<std::size_t>(b + i)]; };
        if (splitAt(0) >= piF) {
            return 0;
        }
        int vote = 0;
        for (int i = 0; i < 4; ++i) {
            const bool flat = angles->flatterUnsplit[static_cast<std::size_t>(topology.vertex(b + i))] != 0;
            vote += (i % 2 == 0) ? (flat ? -1 : 1) : (flat ? 1 : -1);
            if (i == 3) {
                return vote >= 0 ? 0 : 1;
            }
            if (splitAt(i + 1) >= piF) {
                return i + 1;
            }
        }
    }
    const auto& points = mesh.points();
    auto p = [&](int i) { return points[static_cast<std::size_t>(topology.vertex(b + i))]; };
    auto squared = [](Vec3f d) { return (d.x * d.x + d.y * d.y) + d.z * d.z; };
    return squared(p(0) - p(2)) > squared(p(1) - p(3)) ? 1 : 0;
}

// Start corner of a face on the mirror path.
int mirroredStart(const Mesh& mesh, const MirrorLists& mirror, std::vector<int>& startPoint, int f)
{
    const Topology& topology = mesh.topology();
    const int b = topology.faceBegin(f);
    const int n = topology.faceSize(f);
    int start = 0;
    const int m = static_cast<std::size_t>(f) < mirror.faces.size() ? mirror.faces[static_cast<std::size_t>(f)] : -1;
    if (m >= 0 && m < f) {
        const int sp = startPoint[static_cast<std::size_t>(m)];
        const int target = sp >= 0 && static_cast<std::size_t>(sp) < mirror.points.size()
                               ? mirror.points[static_cast<std::size_t>(sp)]
                               : -1;
        for (int i = 1; i < n; ++i) {
            if (topology.vertex(b + i) == target) {
                start = i;
                break;
            }
        }
    } else {
        const int two = twoCornerStart(mesh.adjacency(), topology, f);
        start = two > 0 ? two : 0;
    }
    startPoint[static_cast<std::size_t>(f)] = topology.vertex(b + start);
    return start;
}

} // namespace

Mesh Triangulator::run(Mesh&& source, const Options& options)
{
    if (source.topology().mode() == FaceMode::Triangles) {
        return std::move(source);
    }
    const Topology& topology = source.topology();
    const bool mirrored = options.mirror != nullptr && !options.mirror->faces.empty();
    // With the mirror switch on the angle rule stays off, also when the mesh
    // has no face pairs: quads then take the shorter diagonal.
    std::optional<CornerAngles> angles;
    if (options.angleRule && options.mirror == nullptr) {
        angles = cornerAngles(source);
    }
    const std::optional<BitArray> sharp = SharpEdges::halfEdgeBits(source);
    const auto& points = source.points();

    std::vector<int> corners;      // point per new corner
    std::vector<int> cornerOrigin; // source corner per new corner
    std::vector<int> edgeOrigin;   // source corner whose half edge the new one is, or -1
    std::vector<int> faceOrigin;   // source face per triangle
    std::vector<int> startPoint(mirrored ? static_cast<std::size_t>(topology.faceCount()) : 0, -1);
    LocalTriangles local;
    std::vector<Vec3f> polygon;
    for (int f = 0; f < topology.faceCount(); ++f) {
        const int b = topology.faceBegin(f);
        const int n = topology.faceSize(f);
        local.clear();
        if (mirrored) {
            PolygonTriangulator::fan(n, mirroredStart(source, *options.mirror, startPoint, f), local);
        } else if (n == 4) {
            PolygonTriangulator::fan(4, quadStart(source, angles ? &*angles : nullptr, f), local);
        } else if (n == 3) {
            local.push_back({0, 1, 2});
        } else if (n > 4) {
            polygon.clear();
            for (int i = 0; i < n; ++i) {
                polygon.push_back(points[static_cast<std::size_t>(topology.vertex(b + i))]);
            }
            local = PolygonTriangulator::triangulate(polygon);
        }
        for (const auto& t : local) {
            for (int j = 0; j < 3; ++j) {
                const int i = t[static_cast<std::size_t>(j)];
                const int following = t[static_cast<std::size_t>((j + 1) % 3)];
                corners.push_back(topology.vertex(b + i));
                cornerOrigin.push_back(b + i);
                edgeOrigin.push_back((i + 1) % n == following ? b + i : -1);
            }
            faceOrigin.push_back(f);
        }
    }

    Mesh out;
    std::vector<Vec3f> movedPoints;
    source.swapPoints(movedPoints);
    out.swapPoints(movedPoints);
    out.editTopology().swapCorners(corners, FaceMode::Triangles);
    out.topologyChanged();
    for (const UvSet& set : source.uvSets()) {
        UvSet uv;
        uv.name = set.name;
        uv.coords = set.coords;
        if (!set.cornerUv.empty()) {
            uv.cornerUv.reserve(cornerOrigin.size());
            for (const int c : cornerOrigin) {
                uv.cornerUv.push_back(set.cornerUv[static_cast<std::size_t>(c)]);
            }
        }
        out.uvSets().push_back(std::move(uv));
    }
    IndexMap faces;
    faces.newToOld = std::move(faceOrigin);
    for (std::size_t k = 0; k < source.maps().count(); ++k) {
        const AttributeMap& map = source.maps().at(k);
        if (map.kind() == ElementKind::Point) {
            out.maps().set(map.clone());
        } else if (map.kind() == ElementKind::Face) {
            std::unique_ptr<AttributeMap> copy = map.clone();
            copy->reorder(faces);
            out.maps().set(std::move(copy));
        }
    }
    if (sharp) {
        BitArray bits(edgeOrigin.size());
        if (!mirrored) {
            for (std::size_t c = 0; c < edgeOrigin.size(); ++c) {
                const int o = edgeOrigin[c];
                if (o >= 0 && static_cast<std::size_t>(o) < sharp->size() && sharp->test(static_cast<std::size_t>(o))) {
                    bits.set(c);
                }
            }
        }
        out.maps().set(std::make_unique<BoolMap>(std::string(map_names::sharpEdges), ElementKind::HalfEdge, std::move(bits)));
    }
    source = Mesh();
    return out;
}

} // namespace areko::remesher::prep
