#include "net_lines.h"

#include "face_borders.h"
#include "mesh/edge_chains.h"
#include "mesh/link_path.h"
#include "mesh/sharp_edges.h"

#include <algorithm>
#include <iterator>
#include <numeric>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

bool hasSharpEdges(const Mesh& mesh)
{
    return mesh.maps().find(map_names::sharpEdges, ElementKind::Edge) != nullptr ||
           mesh.maps().find(map_names::sharpEdges, ElementKind::HalfEdge) != nullptr;
}

void setEdge(BitArray* edges, int e)
{
    if (edges != nullptr) {
        edges->setBounded(static_cast<std::size_t>(e));
    }
}

} // namespace

void NetLines::reset(bool hasEdgeAttribute)
{
    state_ = hasEdgeAttribute ? 0 : -1;
    lines_.clear();
    recordedPointCount_ = 0;
}

std::vector<NetLine> NetLines::release()
{
    std::vector<NetLine> out = std::move(lines_);
    reset(false);
    return out;
}

void NetLines::setIdentityOrigins(Mesh& mesh)
{
    std::vector<int> identity(static_cast<std::size_t>(mesh.pointCount()));
    std::iota(identity.begin(), identity.end(), 0);
    mesh.maps().set(std::make_unique<TypedMap<int>>(std::string(originMap), ElementKind::Point, std::move(identity)));
    recordedPointCount_ = mesh.pointCount();
}

void NetLines::extractFromSharpEdges(Mesh& mesh)
{
    if (state_ != 0 || mesh.topology().mode() != FaceMode::Triangles) {
        return;
    }
    if (!hasSharpEdges(mesh) && mesh.maps().faceGroups() == nullptr) {
        return;
    }
    lines_ = linesFromEdges(mesh, collectEdges(mesh, false, true, false), 1);
    state_ = 1;
    setIdentityOrigins(mesh);
    mesh.maps().remove(map_names::sharpEdges, ElementKind::HalfEdge);
    mesh.maps().remove(map_names::sharpEdges, ElementKind::Edge);
}

void NetLines::appendChainLines(const Mesh& mesh, const Csr& chains, const ChainForm& form,
                                std::vector<NetLine>& out)
{
    const auto& points = mesh.points();
    for (int r = 0; r < chains.rows(); ++r) {
        NetLine line;
        for (const int v : chains.row(r)) {
            if (form.validOnly && v < 0) {
                continue;
            }
            line.vertices.push_back(v);
            line.points.push_back(points[static_cast<std::size_t>(v)]);
        }
        if (form.validOnly && line.vertices.size() < 2) {
            continue;
        }
        line.type = form.types.empty() ? 2 : form.types[static_cast<std::size_t>(r)];
        line.flags = form.flags;
        out.push_back(std::move(line));
    }
}

std::vector<NetLine> NetLines::linesFromEdges(const Mesh& mesh, const BitArray& edges, std::uint8_t flags)
{
    std::vector<NetLine> out;
    appendChainLines(mesh, EdgeChains::build(mesh.adjacency(), edges), {.flags = flags, .validOnly = true}, out);
    return out;
}

std::optional<NetLines::FinalLines> NetLines::buildFinal(Mesh& mesh, FinalSources sources)
{
    if (mesh.topology().mode() != FaceMode::Triangles) {
        return std::nullopt;
    }
    const TypedMap<int>* groups = mesh.maps().faceGroups();
    const bool anything = state_ == 1 ? (groups != nullptr || !lines_.empty())
                                      : (state_ == 0 ? hasSharpEdges(mesh) : groups != nullptr);
    if (!anything) {
        return std::nullopt;
    }
    const Adjacency& adjacency = mesh.adjacency();
    const auto edgeCount = static_cast<std::size_t>(adjacency.edgeCount());
    FinalLines out{BitArray(edgeCount), {}};
    std::vector<NetLine> lines;
    BitArray sharp(edgeCount);
    if (sources.sharpEdges) {
        sharp = collectEdges(mesh, false, true, false);
        out.protectedEdges.orWith(sharp);
        lines = linesFromEdges(mesh, sharp, 1);
    }
    if (sources.groupBorders) {
        // Without face groups there are no group borders.
        BitArray borders = groups != nullptr ? FaceBorders::valueBorders(adjacency, groups->values(), false)
                                             : BitArray(edgeCount);
        out.protectedEdges.orWith(borders);
        for (std::size_t e = 0; e < borders.size(); ++e) {
            if (sharp.testBounded(e)) {
                borders.set(e, false);
            }
        }
        std::vector<NetLine> groupLines = linesFromEdges(mesh, borders, 2);
        lines.insert(lines.end(), std::make_move_iterator(groupLines.begin()),
                     std::make_move_iterator(groupLines.end()));
    }
    lines_ = std::move(lines);
    state_ = 1;
    out.linePoints.resize(static_cast<std::size_t>(mesh.pointCount()));
    for (const NetLine& line : lines_) {
        for (const int v : line.vertices) {
            out.linePoints.setBounded(static_cast<std::size_t>(v));
        }
    }
    return out;
}

void NetLines::appendChains(const Mesh& mesh, const Csr& chains, std::span<const int> types, std::uint8_t flags)
{
    appendChainLines(mesh, chains, {.types = types, .flags = flags}, lines_);
}

void NetLines::detectedToProtection()
{
    for (NetLine& line : lines_) {
        if (line.flags == 4) {
            line.flags = 8;
            line.type = 3;
        }
    }
}

void NetLines::restore(Mesh& mesh, Restore mode, BitArray* edges)
{
    const auto* origins = mesh.maps().findAs<TypedMap<int>>(originMap, ElementKind::Point);
    if (origins == nullptr || origins->values().empty()) {
        return;
    }
    const std::vector<int>& map = origins->values();
    const int recorded = recordedPointCount_;
    std::vector<int> current(static_cast<std::size_t>(std::max(recorded, 0)), -1);
    for (int i = static_cast<int>(map.size()) - 1; i >= 0; --i) {
        const int o = map[static_cast<std::size_t>(i)];
        if (o >= 0 && o < recorded) {
            current[static_cast<std::size_t>(o)] = i;
        }
    }
    auto newNumber = [&](int old) { return (old >= 0 && old < recorded) ? current[static_cast<std::size_t>(old)] : -1; };

    const Adjacency& adjacency = mesh.adjacency();
    const auto& points = mesh.points();
    if (edges != nullptr) {
        edges->resize(static_cast<std::size_t>(adjacency.edgeCount()));
    }
    for (NetLine& line : lines_) {
        std::size_t k = 0;
        while (k < line.vertices.size() && newNumber(line.vertices[k]) < 0) {
            ++k;
        }
        if (k == line.vertices.size()) {
            continue;
        }
        std::vector<int> vertices{newNumber(line.vertices[k])};
        int previous = vertices.front();
        for (++k; k < line.vertices.size(); ++k) {
            const int b = newNumber(line.vertices[k]);
            if (b == -1) {
                continue;
            }
            if (mode == Restore::OverTwoEdges) {
                const LinkPath path = shortestLink(adjacency, points, previous, b);
                setEdge(edges, path.first);
                if (path.via >= 0) {
                    setEdge(edges, path.second);
                    vertices.push_back(path.via);
                }
            } else {
                setEdge(edges, adjacency.edgeBetween(previous, b));
            }
            vertices.push_back(b);
            previous = b;
        }
        line.points.clear();
        for (const int v : vertices) {
            line.points.push_back(points[static_cast<std::size_t>(v)]);
        }
        line.vertices = std::move(vertices);
    }
    setIdentityOrigins(mesh);
}

void NetLines::toSharpEdges(Mesh& mesh)
{
    if (state_ != 1) {
        return;
    }
    BitArray bits(static_cast<std::size_t>(mesh.adjacency().edgeCount()));
    addEdges(mesh, lines_, true, bits);
    std::erase_if(lines_, isSharpLine);
    SharpEdges::setEdgeBits(mesh, std::move(bits));
    state_ = 0;
}

void NetLines::addEdges(const Mesh& mesh, std::span<const NetLine> lines, bool onlySharp, BitArray& edges)
{
    const Adjacency& adjacency = mesh.adjacency();
    const int pointCount = mesh.pointCount();
    for (const NetLine& line : lines) {
        if (onlySharp && !isSharpLine(line)) {
            continue;
        }
        for (std::size_t i = 1; i < line.vertices.size(); ++i) {
            const int a = line.vertices[i - 1];
            const int b = line.vertices[i];
            if (a >= 0 && b >= 0 && a < pointCount && b < pointCount) {
                setEdge(&edges, adjacency.edgeBetween(a, b));
            }
        }
    }
}

BitArray NetLines::collectEdges(Mesh& mesh, bool groups, bool sharp, bool withBoundary, Restore restoreMode)
{
    BitArray edges(static_cast<std::size_t>(mesh.adjacency().edgeCount()));
    if (groups) {
        if (const auto* values = mesh.maps().faceGroups()) {
            edges.orWith(FaceBorders::valueBorders(mesh.adjacency(), values->values(), withBoundary));
        } else {
            for (const int e : mesh.adjacency().boundaryEdges()) {
                edges.set(static_cast<std::size_t>(e));
            }
        }
    }
    if (sharp) {
        if (state_ == 1) {
            restore(mesh, restoreMode, &edges);
        } else if (state_ == 0) {
            if (const auto bits = SharpEdges::edgeBits(mesh)) {
                edges.orWith(*bits);
            }
        }
    }
    return edges;
}

} // namespace areko::remesher::prep
