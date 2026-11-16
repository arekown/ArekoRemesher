#include "part_builder.h"

#include "mesh/connected_parts.h"
#include "mesh/index_map.h"

#include <cstddef>
#include <utility>

namespace areko::remesher::parts {

using namespace areko::remesher::mesh;
using prep::SymmetryCache;

namespace {

constexpr double steepAngle = 1.7453292519943295; // 100 degrees

// Point and face maps of `source` for the elements listed new -> old.
AttributeMaps subsetMaps(const AttributeMaps& source, const std::vector<int>& points, const std::vector<int>& faces)
{
    AttributeMaps maps(source);
    for (const ElementKind kind : {ElementKind::Corner, ElementKind::UvPoint, ElementKind::Edge, ElementKind::UvEdge,
                                   ElementKind::HalfEdge}) {
        maps.removeKind(kind);
    }
    IndexMap pointMap;
    pointMap.newToOld = points;
    pointMap.count = static_cast<int>(points.size());
    IndexMap faceMap;
    faceMap.newToOld = faces;
    faceMap.count = static_cast<int>(faces.size());
    for (std::size_t i = 0; i < maps.count(); ++i) {
        AttributeMap& map = maps.at(i);
        (void)map.reorder(map.kind() == ElementKind::Point ? pointMap : faceMap);
    }
    return maps;
}

} // namespace

PartBuilder::Split PartBuilder::split(const Mesh& prepared, const SymmetryCache* mirror, bool symmetryOn)
{
    Split out;
    const ConnectedParts& components = prepared.connectedParts();
    const int pointCount = prepared.pointCount();
    for (int v = 0; v < pointCount; ++v) {
        if (components.partOfPoint(v) == -2) {
            out.sharedPoints = true;
            return out;
        }
    }
    const int partCount = components.partCount();
    const Csr& partFaces = components.partFaces();
    out.faceParts = components.faceParts();
    out.facePlaces.assign(static_cast<std::size_t>(prepared.faceCount()), -1);

    // Local number of every point: running count within its component.
    std::vector<int> local(static_cast<std::size_t>(pointCount), -1);
    std::vector<std::vector<int>> partPoints(static_cast<std::size_t>(partCount));
    for (int v = 0; v < pointCount; ++v) {
        const int p = components.partOfPoint(v);
        if (p < 0) {
            continue;
        }
        auto& list = partPoints[static_cast<std::size_t>(p)];
        local[static_cast<std::size_t>(v)] = static_cast<int>(list.size());
        list.push_back(v);
    }

    const Topology& topology = prepared.topology();
    const bool withMirror = symmetryOn && mirror != nullptr && !mirror->pointMirror().empty();
    out.parts.resize(static_cast<std::size_t>(partCount));
    for (int p = 0; p < partCount; ++p) {
        Part& part = out.parts[static_cast<std::size_t>(p)];
        part.index = p;
        const std::vector<int>& points = partPoints[static_cast<std::size_t>(p)];
        const auto faces = partFaces.row(p);
        std::vector<Vec3f> positions;
        positions.reserve(points.size());
        for (const int v : points) {
            positions.push_back(prepared.points()[static_cast<std::size_t>(v)]);
        }
        std::vector<int> corners;
        corners.reserve(3 * faces.size());
        std::vector<int> faceList(faces.begin(), faces.end());
        for (std::size_t i = 0; i < faces.size(); ++i) {
            const int f = faces[i];
            out.facePlaces[static_cast<std::size_t>(f)] = static_cast<int>(i);
            for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
                corners.push_back(local[static_cast<std::size_t>(topology.vertex(c))]);
            }
        }
        Mesh& mesh = part.mesh;
        mesh.swapPoints(positions);
        mesh.editTopology().swapCorners(corners, FaceMode::Triangles);
        mesh.topologyChanged();
        mesh.maps() = subsetMaps(prepared.maps(), points, faceList);
        mesh.setReplacementNormalMode(2, true);

        if (withMirror) {
            const std::vector<int>& pairs = mirror->pointMirror();
            const BitArray& sides = mirror->pointSide();
            std::vector<int> localPairs(points.size(), -1);
            BitArray localSides(sides.empty() ? 0 : points.size());
            for (std::size_t i = 0; i < points.size(); ++i) {
                const int v = points[i];
                const int m = static_cast<std::size_t>(v) < pairs.size() ? pairs[static_cast<std::size_t>(v)] : -1;
                if (m >= 0 && m < pointCount && components.partOfPoint(m) == p) {
                    localPairs[i] = local[static_cast<std::size_t>(m)];
                }
                if (!localSides.empty()) {
                    localSides.set(i, sides.testBounded(static_cast<std::size_t>(v)));
                }
            }
            part.mirror = PartMirror(SymmetryCache::fromPoints(std::move(localPairs), std::move(localSides),
                                                               mirror->lockedPlane()));
        }
    }
    return out;
}

void PartBuilder::assignBudgets(std::vector<Part>& parts, const Mesh& prepared, int quadBudget, float edgeLength)
{
    if (parts.empty()) {
        return;
    }
    const float total = prepared.triangleGeometry().totalArea();
    const float n = !(edgeLength > 0.0f) ? static_cast<float>(quadBudget) : total / (edgeLength * edgeLength);
    if (parts.size() == 1) {
        parts.front().quadBudget = truncToInt(n);
        parts.front().targetEdgeLength = edgeLength;
        return;
    }
    const float perArea = n / total;
    for (Part& part : parts) {
        const float x = part.mesh.triangleGeometry().totalArea() * perArea;
        part.quadBudget = x < 10.0f ? 10 : truncToInt(x);
        part.targetEdgeLength = edgeLength;
    }
}

bool PartBuilder::negligible(const Part& part, const Mesh& prepared)
{
    const Mesh& mesh = part.mesh;
    const int points = mesh.pointCount();
    if (points == 3 && mesh.cornerCount() == 6) {
        return true;
    }
    const float share = mesh.triangleGeometry().totalArea() / prepared.triangleGeometry().totalArea();
    const double r = static_cast<double>(share);
    if (r < 2e-6) {
        if (points < 7) {
            return true;
        }
        if (points > 20) {
            return false;
        }
        return r < 1e-9;
    }
    if (!(r < 1e-4) || points > 5) {
        return false;
    }
    const Adjacency& adjacency = mesh.adjacency();
    const std::vector<float>& dihedral = mesh.dihedralAngles();
    const int edges = adjacency.edgeCount();
    int steep = 0;
    int border = 0;
    for (int e = 0; e < edges; ++e) {
        if (static_cast<double>(dihedral[static_cast<std::size_t>(e)]) > steepAngle) {
            ++steep;
        }
        if (adjacency.isBoundaryEdge(e)) {
            ++border;
        }
    }
    if (border > 0) {
        return edges == steep + border;
    }
    return steep > edges / 2 || edges == steep;
}

const Mesh* PartBuilder::rawMeshOf(const Part& part)
{
    return part.raw ? &part.raw->mesh : nullptr;
}

Mesh PartBuilder::merged(const std::vector<Part>& parts, MeshOf meshOf, std::vector<int>& faceParts,
                         std::vector<int>& pointStarts, MapRules rules)
{
    Mesh out;
    std::vector<Vec3f> points;
    std::vector<int> faceStart{0};
    std::vector<int> corners;
    faceParts.clear();
    pointStarts.assign(1, 0);
    bool first = true;
    bool triangles = true;
    int edges = 0;
    for (const Part& part : parts) {
        if (const Mesh* mesh = meshOf(part)) {
            const ElementCounts own{static_cast<int>(points.size()), static_cast<int>(faceStart.size()) - 1,
                                    static_cast<int>(corners.size()), edges};
            const bool empty = rules.result ? own.points == 0 && own.faces == 0 : first;
            const int offset = static_cast<int>(points.size());
            points.insert(points.end(), mesh->points().begin(), mesh->points().end());
            const Topology& topology = mesh->topology();
            triangles = triangles && topology.mode() == FaceMode::Triangles;
            for (int f = 0; f < topology.faceCount(); ++f) {
                for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
                    corners.push_back(topology.vertex(c) + offset);
                }
                faceStart.push_back(static_cast<int>(corners.size()));
            }
            faceParts.insert(faceParts.end(), static_cast<std::size_t>(topology.faceCount()), part.index);
            const int meshEdges = rules.result ? mesh->adjacency().edgeCount() : 0;
            if (empty) {
                out.maps() = mesh->maps();
                first = false;
            } else {
                out.maps().appendFrom(mesh->maps(),
                                      ElementCounts{mesh->pointCount(), mesh->faceCount(), mesh->cornerCount(), meshEdges},
                                      rules.unite ? &own : nullptr);
            }
            edges += meshEdges;
        }
        pointStarts.push_back(static_cast<int>(points.size()));
    }
    out.swapPoints(points);
    if (triangles) {
        out.editTopology().swapCorners(corners, FaceMode::Triangles);
    } else {
        out.editTopology().swapPolygons(faceStart, corners);
    }
    out.topologyChanged();
    return out;
}

} // namespace areko::remesher::parts
