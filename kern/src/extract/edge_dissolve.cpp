#include "edge_dissolve.h"

#include "mesh/bit_array.h"
#include "mesh/connected_parts.h"
#include "mesh/edge_chains.h"
#include "mesh/measures.h"
#include "mesh/modification.h"

#include <algorithm>
#include <vector>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

// Groups of faces and the points left out, built from one call.
class Dissolution {
public:
    Dissolution(const Mesh& mesh, std::span<const int> edges)
        : mesh_(mesh), adjacency_(mesh.adjacency()), topology_(mesh.topology())
    {
        const int edgeCount = adjacency_.edgeCount();
        const int faceCount = topology_.faceCount();
        affected_.resize(static_cast<std::size_t>(faceCount));
        for (const int e : edges) {
            for (const int c : adjacency_.edgeCorners(e)) {
                affected_.set(static_cast<std::size_t>(adjacency_.cornerFace(c)));
            }
        }
        // Faces of the list joined over the listed edges only; a face across
        // a listed edge is always a listed face.
        BitArray blocked(static_cast<std::size_t>(edgeCount), true);
        for (const int e : edges) {
            blocked.set(static_cast<std::size_t>(e), false);
        }
        int components = 0;
        const std::vector<int> component = ConnectedParts::faceComponents(topology_, adjacency_, &blocked, components);
        std::vector<int> number(static_cast<std::size_t>(components), -1);
        group_.assign(static_cast<std::size_t>(faceCount), -1);
        for (int f = 0; f < faceCount; ++f) {
            if (!affected_.test(static_cast<std::size_t>(f))) {
                continue;
            }
            int& g = number[static_cast<std::size_t>(component[static_cast<std::size_t>(f)])];
            if (g < 0) {
                g = static_cast<int>(groupFaces_.size());
                groupFaces_.emplace_back();
            }
            group_[static_cast<std::size_t>(f)] = g;
            groupFaces_[static_cast<std::size_t>(g)].push_back(f);
        }
        removed_.resize(static_cast<std::size_t>(mesh.pointCount()));
    }

    // Leaves out the points whose two remaining edges are nearly straight;
    // their faces outside the groups become groups of their own.
    void removeStraightPoints()
    {
        const int edgeCount = adjacency_.edgeCount();
        std::vector<char> inner(static_cast<std::size_t>(edgeCount), 0);
        for (int e = 0; e < edgeCount; ++e) {
            const auto corners = adjacency_.edgeCorners(e);
            if (corners.size() < 2) {
                continue;
            }
            const int g = group_[static_cast<std::size_t>(adjacency_.cornerFace(corners[0]))];
            bool same = g >= 0;
            for (std::size_t k = 1; k < corners.size() && same; ++k) {
                same = group_[static_cast<std::size_t>(adjacency_.cornerFace(corners[k]))] == g;
            }
            inner[static_cast<std::size_t>(e)] = same ? 1 : 0;
        }

        BitArray candidates(static_cast<std::size_t>(mesh_.pointCount()));
        affected_.forEachSet([&](std::size_t f) {
            for (int c = topology_.faceBegin(static_cast<int>(f)); c < topology_.faceEnd(static_cast<int>(f)); ++c) {
                if (const int v = topology_.vertex(c); v >= 0) {
                    candidates.set(static_cast<std::size_t>(v));
                }
            }
        });

        const auto& points = mesh_.points();
        std::vector<int> extra;
        candidates.forEachSet([&](std::size_t p) {
            const int v = static_cast<int>(p);
            int rest[2] = {-1, -1};
            int count = 0;
            for (const int e : adjacency_.pointEdges(v)) {
                if (inner[static_cast<std::size_t>(e)]) {
                    continue;
                }
                if (count < 2) {
                    rest[count] = e;
                }
                if (++count == 3) {
                    break;
                }
            }
            if (count != 2) {
                return;
            }
            const int a = adjacency_.otherEnd(rest[0], v);
            const int b = adjacency_.otherEnd(rest[1], v);
            if (a < 0 || b < 0) {
                return;
            }
            const Vec3f at = points[p];
            const float bend = vectorAngle(at - points[static_cast<std::size_t>(a)], points[static_cast<std::size_t>(b)] - at);
            if (!(bend <= EdgeDissolve::bendLimit)) {
                return;
            }
            removed_.set(p);
            for (const int c : adjacency_.pointCorners(v)) {
                const int f = adjacency_.cornerFace(c);
                if (!affected_.test(static_cast<std::size_t>(f))) {
                    affected_.set(static_cast<std::size_t>(f));
                    extra.push_back(f);
                }
            }
        });
        for (const int f : extra) {
            group_[static_cast<std::size_t>(f)] = static_cast<int>(groupFaces_.size());
            groupFaces_.push_back({f});
        }
    }

    // One polygon per border chain of every group, in group order.
    void buildPolygons(std::vector<int>& faceStart, std::vector<int>& corners, std::vector<int>& faceSource) const
    {
        const int edgeCount = adjacency_.edgeCount();
        BitArray inGroup(static_cast<std::size_t>(topology_.faceCount()));
        std::vector<int> chain;
        for (const std::vector<int>& faces : groupFaces_) {
            for (const int f : faces) {
                inGroup.set(static_cast<std::size_t>(f));
            }
            BitArray border(static_cast<std::size_t>(edgeCount));
            for (const int f : faces) {
                for (int c = topology_.faceBegin(f); c < topology_.faceEnd(f); ++c) {
                    const int e = adjacency_.cornerEdge(c);
                    if (e >= 0 && isGroupBorder(e, inGroup)) {
                        border.set(static_cast<std::size_t>(e));
                    }
                }
            }
            const Csr chains = EdgeChains::build(adjacency_, border);
            for (int r = 0; r < chains.rows(); ++r) {
                const auto row = chains.row(r);
                chain.assign(row.begin(), row.end());
                if (chain.size() >= 2) {
                    align(chain, inGroup);
                }
                if (chain.front() == chain.back()) {
                    chain.pop_back();
                }
                const std::size_t before = corners.size();
                for (const int v : chain) {
                    if (!removed_.testBounded(static_cast<std::size_t>(v))) {
                        corners.push_back(v);
                    }
                }
                if (corners.size() - before >= 3) {
                    faceStart.push_back(static_cast<int>(corners.size()));
                    faceSource.push_back(faces.front());
                } else {
                    corners.resize(before);
                }
            }
            for (const int f : faces) {
                inGroup.set(static_cast<std::size_t>(f), false);
            }
        }
    }

    [[nodiscard]] const BitArray& affected() const noexcept { return affected_; }

private:
    // An edge with two corners and one face in the group, a border edge of a
    // face in the group, or an edge with more corners and faces in and outside.
    [[nodiscard]] bool isGroupBorder(int e, const BitArray& inGroup) const
    {
        const auto corners = adjacency_.edgeCorners(e);
        auto in = [&](int c) { return inGroup.test(static_cast<std::size_t>(adjacency_.cornerFace(c))); };
        if (corners.size() == 2) {
            return in(corners[0]) != in(corners[1]);
        }
        if (corners.size() == 1) {
            return in(corners[0]);
        }
        bool inside = false;
        bool outside = false;
        for (const int c : corners) {
            (in(c) ? inside : outside) = true;
        }
        return inside && outside;
    }

    // Turns the chain round when the group's first face at its first edge
    // runs the other way.
    void align(std::vector<int>& chain, const BitArray& inGroup) const
    {
        const int e = adjacency_.edgeBetween(chain[0], chain[1]);
        if (e < 0) {
            return;
        }
        for (const int c : adjacency_.edgeCorners(e)) {
            if (!inGroup.test(static_cast<std::size_t>(adjacency_.cornerFace(c)))) {
                continue;
            }
            if (adjacency_.vertex(c) == chain[1] && adjacency_.vertex(adjacency_.next(c)) == chain[0]) {
                std::reverse(chain.begin(), chain.end());
            }
            return;
        }
    }

    const Mesh& mesh_;
    const Adjacency& adjacency_;
    const Topology& topology_;
    BitArray affected_;
    std::vector<int> group_;
    std::vector<std::vector<int>> groupFaces_;
    BitArray removed_;
};

} // namespace

void EdgeDissolve::run(Mesh& mesh, std::span<const int> edges, bool removeStraightPoints)
{
    mesh.toPolygonMode();
    std::vector<int> faceStart{0};
    std::vector<int> corners;
    std::vector<int> faceSource;
    BitArray affected;
    {
        Dissolution dissolution(mesh, edges);
        if (removeStraightPoints) {
            dissolution.removeStraightPoints();
        }
        dissolution.buildPolygons(faceStart, corners, faceSource);
        affected = dissolution.affected();
    }
    // The new polygons refer to the points of the mesh and bring none of their own.
    Mesh appended;
    appended.editTopology().swapPolygons(faceStart, corners);
    Modification change;
    change.deleteFaces = std::move(affected);
    change.appendFaces = &appended;
    if (mesh.maps().hasKind(ElementKind::Face)) {
        change.faces.newToOld = std::move(faceSource);
    }
    MeshModifier::apply(mesh, change);
}

} // namespace areko::remesher::extract
