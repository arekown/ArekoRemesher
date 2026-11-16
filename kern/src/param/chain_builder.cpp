#include "chain_builder.h"

#include "chain_sides.h"
#include "curve_chain_builder.h"

#include "mesh/edge_chains.h"

#include <algorithm>
#include <array>
#include <utility>

namespace areko::remesher::param {

using namespace areko::remesher::mesh;

namespace {

// The working state of one computation.
class Computation {
public:
    Computation(parts::Part& part, const ChainBuilder::Options& options, const run::Progress& progress)
        : part_(part), mesh_(part.mesh), adjacency_(part.mesh.adjacency()), graph_(*part.cutGraph), options_(options),
          progress_(progress)
    {
    }

    [[nodiscard]] Chains run();

private:
    void chainCuts();
    // The node that stands in for the missing singular points, or -1.
    [[nodiscard]] int substituteNode();
    void markSingularNodes();
    void markFullTurnChains();
    void followSides();
    void markCurveCheckChains();
    void assignOffsetVariables();

    parts::Part& part_;
    const Mesh& mesh_;
    const Adjacency& adjacency_;
    const field::CutGraph& graph_;
    const ChainBuilder::Options& options_;
    const run::Progress& progress_;
    Chains out_;
};

Chains Computation::run()
{
    out_.singularities = part_.crossField->singularities;
    chainCuts();
    if (out_.singularities.points.empty() && !out_.nodes.empty()) {
        if (const int node = substituteNode(); node != -1) {
            out_.singularities.points.push_back(node);
            out_.singularities.index.push_back(0);
            out_.singularities.fullTurn.push_back(0);
            out_.substituted = true;
        }
    }
    markSingularNodes();
    markFullTurnChains();
    progress_.check();
    followSides();
    progress_.check();
    markCurveCheckChains();
    assignOffsetVariables();
    if (!graph_.curveEdges.empty()) {
        out_.curves = CurveChainBuilder::build(mesh_, *part_.crossField, graph_, options_.symmetry);
    }
    progress_.check();
    return std::move(out_);
}

void Computation::chainCuts()
{
    EdgeChains::Options rules;
    rules.separators = out_.singularities.points;
    rules.borderPointsSeparate = true;
    rules.startsFollowChains = true;
    EdgeChains::Result cuts = EdgeChains::build(adjacency_, graph_.cutEdgeBits, rules);
    out_.points = std::move(cuts.points);
    out_.edges = std::move(cuts.edges);
    out_.nodes = std::move(cuts.nodes);
}

int Computation::substituteNode()
{
    if (!options_.symmetry) {
        return out_.nodes.front();
    }
    if (part_.mirror.symmetric(mesh_, true) == nullptr || graph_.mirrorHalfCount == 0) {
        return -1;
    }
    // The first node with a triangle in the computed half.
    for (const int node : out_.nodes) {
        const auto corners = adjacency_.pointCorners(node);
        if (std::any_of(corners.begin(), corners.end(),
                        [&](int c) { return !graph_.mirrorHalf.testBounded(static_cast<std::size_t>(c / 3)); })) {
            return node;
        }
    }
    return -1;
}

void Computation::markSingularNodes()
{
    const std::vector<int>& singular = out_.singularities.points;
    out_.singularNodes = BitArray(out_.nodes.size());
    for (std::size_t i = 0; i < out_.nodes.size(); ++i) {
        if (std::find(singular.begin(), singular.end(), out_.nodes[i]) != singular.end()) {
            out_.singularNodes.set(i);
        }
    }
}

void Computation::markFullTurnChains()
{
    out_.fullTurnChains = BitArray(static_cast<std::size_t>(out_.count()));
    const field::Singularities& singular = out_.singularities;
    std::vector<int> wholeTurns;
    for (std::size_t i = 0; i < singular.points.size(); ++i) {
        if ((singular.fullTurn[i] & 1u) != 0) {
            wholeTurns.push_back(singular.points[i]);
        }
    }
    if (wholeTurns.empty()) {
        return;
    }
    auto wholeTurn = [&](int p) { return std::find(wholeTurns.begin(), wholeTurns.end(), p) != wholeTurns.end(); };
    for (int c = 0; c < out_.count(); ++c) {
        if (out_.edges.rowSize(c) != 1) {
            continue;
        }
        const std::array<int, 2> ends = adjacency_.edgeVertices(out_.edges.row(c)[0]);
        int other;
        if (wholeTurn(ends[0])) {
            other = ends[1];
        } else if (wholeTurn(ends[1])) {
            other = ends[0];
        } else {
            continue;
        }
        if (graph_.copy.pointCopies(adjacency_, other).size() == 1) {
            out_.fullTurnChains.set(static_cast<std::size_t>(c));
        }
    }
}

void Computation::followSides()
{
    const std::size_t count = static_cast<std::size_t>(out_.count());
    out_.sideA.assign(count, std::vector<int>{});
    out_.sideB.assign(count, std::vector<int>{});
    out_.rotation.assign(count, 0.0f);
    for (std::size_t c = 0; c < count; ++c) {
        const int chain = static_cast<int>(c);
        ChainSides::Sides sides = ChainSides::of(adjacency_, graph_.copy, graph_.residualJump, out_.edges.row(chain),
                                                 out_.points.row(chain));
        if (sides.found) {
            out_.sideA[c] = std::move(sides.a);
            out_.sideB[c] = std::move(sides.b);
            out_.rotation[c] = static_cast<float>(sides.jump);
        }
    }
}

void Computation::markCurveCheckChains()
{
    out_.curveCheckChains = BitArray(static_cast<std::size_t>(out_.count()));
    for (int c = 0; c < out_.count(); ++c) {
        if (out_.edges.rowSize(c) != 1) {
            continue;
        }
        const int turns = static_cast<int>(out_.rotation[static_cast<std::size_t>(c)]);
        const int edge = out_.edges.row(c)[0];
        if ((turns + 400) % 4 == 0 && graph_.curveCheckCutBits.testBounded(static_cast<std::size_t>(edge))) {
            out_.curveCheckChains.set(static_cast<std::size_t>(c));
        }
    }
}

void Computation::assignOffsetVariables()
{
    out_.offsetVariable.assign(static_cast<std::size_t>(out_.count()), 0);
    // The offsets follow the two coordinates of every point of the cut copy.
    int next = 2 * graph_.copy.pointCount();
    for (int c = 0; c < out_.count(); ++c) {
        if (out_.edges.rowSize(c) == 0 ||
            !graph_.diskCutBits.testBounded(static_cast<std::size_t>(out_.edges.row(c)[0]))) {
            continue;
        }
        out_.genusChains.push_back(c);
        out_.offsetVariable[static_cast<std::size_t>(c)] = next;
        next += 2;
    }
}

} // namespace

void ChainBuilder::compute(parts::Part& part, const Options& options, const run::Progress& progress)
{
    part.chains = Computation(part, options, progress).run();
}

} // namespace areko::remesher::param
