#include "curve_chain_builder.h"

#include "field/curve_alignment.h"
#include "mesh/edge_chains.h"
#include "mesh/vertex_fan.h"

#include <array>
#include <span>

namespace areko::remesher::param {

using namespace areko::remesher::mesh;

namespace {

// The working state of one build.
class Build {
public:
    Build(const Mesh& mesh, const field::CrossField& field, const field::CutGraph& graph, bool symmetry)
        : adjacency_(mesh.adjacency()), copy_(graph.copy.adjacency()), graph_(graph), symmetry_(symmetry),
          labels_(mesh, field, graph), excluded_(static_cast<std::size_t>(mesh.adjacency().edgeCount()))
    {
    }

    [[nodiscard]] CurveChains run();

private:
    void exclude();
    void chainRemaining();
    // Splits the runs `runEdges` (their points in `runPoints`) by label.
    void split(const Csr& runEdges, const Csr& runPoints);
    void append(std::span<const int> edges, int label);
    // Whether a chain with this label runs along the first combed direction.
    [[nodiscard]] static bool alongFirst(int label) noexcept { return label % 2 == 0; }

    const Adjacency& adjacency_;
    const Adjacency& copy_;
    const field::CutGraph& graph_;
    const bool symmetry_;
    field::CurveLabels labels_;
    // Curve edges of the part mesh that get no chain (bit per edge).
    BitArray excluded_;
    CurveChains out_;
    std::vector<int> chainLabels_; // direction label per chain of out_
};

CurveChains Build::run()
{
    exclude();
    chainRemaining();
    return std::move(out_);
}

void Build::exclude()
{
    // A curve edge that runs the same way as the one before it gets no chain.
    const field::CurveLabels::WalkRules rules{
        .startAtOtherBorderEdge = false, .endBeforeReference = false, .checkCutsSeparate = true};
    labels_.curvePoints().forEachSet([&](std::size_t v) {
        labels_.walkAround(static_cast<int>(v), rules, [&](const VertexFan& fan, std::size_t, std::size_t to) {
            excluded_.set(static_cast<std::size_t>(fan.edges[to]));
        });
    });
}

void Build::chainRemaining()
{
    out_.edgeKind.assign(static_cast<std::size_t>(copy_.edgeCount()), field::curve_edge::none);
    // Copy edges along the second combed direction, and along the first.
    std::array<std::vector<int>, 2> lists;
    for (const int e : graph_.curveEdges) {
        if (excluded_.testBounded(static_cast<std::size_t>(e))) {
            continue;
        }
        const std::uint8_t kind = graph_.curveEdgeKind[static_cast<std::size_t>(e)];
        for (const int u : graph_.copy.copyEdges(adjacency_, e)) {
            out_.edgeKind[static_cast<std::size_t>(u)] = kind;
            const int corner = copy_.edgeCorners(u)[0];
            lists[graph_.curveFollowsU.testBounded(static_cast<std::size_t>(corner)) ? 1 : 0].push_back(u);
        }
    }
    for (const std::vector<int>& list : lists) {
        const EdgeChains::Result runs = EdgeChains::build(copy_, list, EdgeChains::Options{});
        split(runs.edges, runs.points);
    }
}

void Build::append(std::span<const int> edges, int label)
{
    out_.edges.entries.insert(out_.edges.entries.end(), edges.begin(), edges.end());
    out_.edges.starts.push_back(static_cast<int>(out_.edges.entries.size()));
    out_.coordinate.push_back(static_cast<std::int8_t>(alongFirst(label) ? 1 : 0));
    chainLabels_.push_back(label);
}

void Build::split(const Csr& runEdges, const Csr& runPoints)
{
    const std::vector<int>& origin = graph_.copy.pointOrigins();
    for (int r = 0; r < runEdges.rows(); ++r) {
        const std::span<const int> edges = runEdges.row(r);
        if (edges.empty()) {
            continue;
        }
        const int firstCorner = copy_.edgeCorners(edges[0])[0];
        if (symmetry_ && graph_.mirrorHalf.testBounded(static_cast<std::size_t>(firstCorner / 3))) {
            continue;
        }
        const std::size_t firstChain = chainLabels_.size();
        std::size_t start = 0;
        int label = -1;
        int w = runPoints.row(r)[0];
        for (std::size_t k = 0; k < edges.size(); ++k) {
            const int next = copy_.otherEnd(edges[k], w);
            const int now = labels_.label(copy_.edgeCorners(edges[k])[0], origin[static_cast<std::size_t>(w)],
                                          origin[static_cast<std::size_t>(next)]);
            if (now != label) {
                if (label != -1) {
                    out_.turnPoints.push_back(origin[static_cast<std::size_t>(w)]);
                    append(edges.subspan(start, k - start), label);
                    start = k;
                }
                label = now;
            }
            w = next;
        }
        const std::span<const int> rest = edges.subspan(start);
        const std::array<int, 2> head = copy_.edgeVertices(edges.front());
        const std::array<int, 2> tail = copy_.edgeVertices(edges.back());
        const bool closed = head[0] == tail[1] || head[0] == tail[0] || tail[1] == head[1] || tail[0] == head[1];
        if (closed && firstChain < chainLabels_.size() && chainLabels_[firstChain] == label) {
            // The run ends where it began: its last piece goes in front of its first.
            std::vector<int>& starts = out_.edges.starts;
            out_.edges.entries.insert(out_.edges.entries.begin() + starts[firstChain], rest.begin(), rest.end());
            for (std::size_t i = firstChain + 1; i < starts.size(); ++i) {
                starts[i] += static_cast<int>(rest.size());
            }
        } else {
            append(rest, label);
        }
    }
}

} // namespace

CurveChains CurveChainBuilder::build(const Mesh& mesh, const field::CrossField& field, const field::CutGraph& graph,
                                     bool symmetry)
{
    return Build(mesh, field, graph, symmetry).run();
}

} // namespace areko::remesher::param
