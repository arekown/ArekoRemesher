#include "crease_split.h"

#include "mesh/sharp_edges.h"
#include "prep/mesh_repair.h"

#include <cstddef>
#include <cstdint>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>

namespace areko::remesher::merge {

using namespace areko::remesher::mesh;

namespace {

class Splitter {
public:
    Splitter(const Mesh& mesh, const BitArray& sharp)
        : adjacency_(mesh.adjacency()), sharp_(sharp), corners_(mesh.topology().cornerVertex()),
          faceLocal_(static_cast<std::size_t>(mesh.faceCount()), 0), next_(mesh.pointCount())
    {
    }

    // Points with at least two cut edges, in the order their second cut edge is counted.
    [[nodiscard]] std::vector<int> candidates() const
    {
        std::vector<std::uint8_t> count(static_cast<std::size_t>(adjacency_.pointCount()), 0);
        std::vector<int> out;
        const auto cut = [&](int e) {
            for (const int v : adjacency_.edgeVertices(e)) {
                std::uint8_t& n = count[static_cast<std::size_t>(v)];
                if (n == 1) {
                    out.push_back(v);
                }
                ++n;
            }
        };
        for (int e = 0; e < adjacency_.edgeCount(); ++e) {
            if (sharp_.testBounded(static_cast<std::size_t>(e))) {
                cut(e);
            }
        }
        for (const int e : adjacency_.boundaryEdges()) {
            if (!sharp_.testBounded(static_cast<std::size_t>(e))) {
                cut(e);
            }
        }
        return out;
    }

    // Gives every fan of v after the first one a new point.
    void split(int v)
    {
        const auto vertexCorners = adjacency_.pointCorners(v);
        const std::size_t n = vertexCorners.size();
        for (std::size_t i = 0; i < n; ++i) {
            faceLocal_[static_cast<std::size_t>(adjacency_.cornerFace(vertexCorners[i]))] = static_cast<int>(i);
        }
        // Cyclic lists of the corners of v (by their face), joined over the edges that are not sharp.
        std::vector<int> next(n);
        std::iota(next.begin(), next.end(), 0);
        std::vector<int> prev = next;
        for (const int e : adjacency_.pointEdges(v)) {
            const auto edgeCorners = adjacency_.edgeCorners(e);
            if (sharp_.testBounded(static_cast<std::size_t>(e)) || edgeCorners.size() < 2) {
                continue;
            }
            const int a = local(edgeCorners[0]);
            int p = prev[static_cast<std::size_t>(a)];
            for (std::size_t k = 1; k < edgeCorners.size(); ++k) {
                const int j = local(edgeCorners[k]);
                const int q = prev[static_cast<std::size_t>(j)];
                next[static_cast<std::size_t>(p)] = j;
                prev[static_cast<std::size_t>(j)] = p;
                next[static_cast<std::size_t>(q)] = a;
                prev[static_cast<std::size_t>(a)] = q;
                p = q;
            }
        }
        // Fans in the order of their smallest corner; the first one keeps v.
        BitArray seen(n);
        bool firstFan = true;
        for (std::size_t i = 0; i < n; ++i) {
            if (seen.test(i)) {
                continue;
            }
            std::size_t j = i;
            do {
                seen.set(j);
                if (!firstFan) {
                    corners_[static_cast<std::size_t>(vertexCorners[j])] = next_;
                }
                j = static_cast<std::size_t>(next[j]);
            } while (j != i);
            if (!firstFan) {
                sources_.push_back(v);
                ++next_;
            }
            firstFan = false;
        }
    }

    std::vector<int>& corners() noexcept { return corners_; }
    [[nodiscard]] const std::vector<int>& sources() const noexcept { return sources_; }

private:
    [[nodiscard]] int local(int corner) const
    {
        return faceLocal_[static_cast<std::size_t>(adjacency_.cornerFace(corner))];
    }

    const Adjacency& adjacency_;
    const BitArray& sharp_;
    std::vector<int> corners_;
    std::vector<int> faceLocal_;
    int next_;
    std::vector<int> sources_;
};

} // namespace

int CreaseSplit::run(Mesh& mesh)
{
    const std::optional<BitArray> sharp = SharpEdges::edgeBits(mesh);
    if (!sharp || !sharp->any()) {
        return 0;
    }
    Splitter splitter(mesh, *sharp);
    for (const int v : splitter.candidates()) {
        splitter.split(v);
    }
    // Applied also without new points.
    prep::MeshRepair::splitPoints(mesh, std::move(splitter.corners()), splitter.sources());
    return static_cast<int>(splitter.sources().size());
}

} // namespace areko::remesher::merge
