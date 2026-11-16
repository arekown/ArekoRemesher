#include "cut_graph_builder.h"

#include "curve_alignment.h"
#include "curve_edges.h"
#include "disk_cut.h"
#include "field_combing.h"
#include "steiner_tree.h"

#include "mesh/boundary_loops.h"
#include "mesh/edge_groups.h"
#include "mesh/measures.h"
#include "prep/freeze_marks.h"
#include "prep/part_symmetry.h"

#include <algorithm>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

namespace {

// The working state of one cut graph computation.
class Computation {
public:
    Computation(parts::Part& part, const CutGraphBuilder::Options& options, const run::Progress& progress)
        : part_(part), mesh_(part.mesh), adjacency_(part.mesh.adjacency()), field_(*part.crossField),
          options_(options), progress_(progress),
          edgeCount_(static_cast<std::size_t>(part.mesh.adjacency().edgeCount()))
    {
    }

    void run();

private:
    void evaluateSymmetry();
    // One pass: cut graph, combing, residual jumps.
    void pass();
    void buildCuts();
    [[nodiscard]] BitArray treeEdges(std::vector<int>& planeChainEdges);
    void addStars();
    void addDiskCuts(const BitArray& rest);
    void mirrorCuts();
    void listCuts();

    parts::Part& part_;
    const Mesh& mesh_;
    const Adjacency& adjacency_;
    const CrossField& field_;
    CutGraphBuilder::Options options_;
    const run::Progress& progress_;
    std::size_t edgeCount_;
    // Mirror pairs of a symmetric part, else null.
    prep::SymmetryCache* symmetric_ = nullptr;
    CutGraph graph_;
};

void Computation::run()
{
    part_.cutGraph.reset();
    evaluateSymmetry();
    progress_.check();
    CurveEdges::collect(part_, part_.mirror.classes.valid ? symmetric_ : nullptr, {options_.followBorders}, graph_);
    pass();
    CurveAlignment::orient(mesh_, field_, graph_);
    CurveAlignment::check(mesh_, field_, graph_);
    if (!graph_.curveCheckCuts.empty()) {
        // The curve check cuts and sources change the tree: once more.
        graph_.rotation.clear();
        graph_.combedU.clear();
        graph_.combedV.clear();
        graph_.residualJump.clear();
        graph_.curveFollowsU.clear();
        pass();
        CurveAlignment::orient(mesh_, field_, graph_);
    }
    progress_.check();
    part_.cutGraph = std::move(graph_);
}

void Computation::evaluateSymmetry()
{
    symmetric_ = part_.mirror.evaluate(mesh_, options_.symmetry);
    progress_.check();
    graph_.mirrorHalf = symmetric_ != nullptr ? prep::PartSymmetry::half(mesh_, *symmetric_, part_.mirror.classes)
                                              : BitArray(static_cast<std::size_t>(mesh_.faceCount()));
    graph_.mirrorHalfCount = static_cast<int>(graph_.mirrorHalf.count());
}

void Computation::pass()
{
    buildCuts();
    progress_.check();
    FieldCombing::comb(mesh_, field_, symmetric_ != nullptr ? &symmetric_->faceMirror(mesh_) : nullptr, graph_);
    progress_.check();
    FieldCombing::residualJumps(adjacency_, field_, graph_);
    progress_.check();
}

void Computation::listCuts()
{
    graph_.cutEdges.clear();
    graph_.cutEdgeBits.forEachSet([&](std::size_t e) { graph_.cutEdges.push_back(static_cast<int>(e)); });
}

void Computation::buildCuts()
{
    // Every pass finds its own disc cuts.
    graph_.diskCuts.clear();
    graph_.diskCutBits.clear();

    std::vector<int> planeChainEdges;
    graph_.cutEdgeBits = treeEdges(planeChainEdges);
    for (const int e : planeChainEdges) {
        graph_.cutEdgeBits.set(static_cast<std::size_t>(e));
    }
    addStars();
    for (const int e : graph_.curveCheckCuts) {
        graph_.cutEdgeBits.set(static_cast<std::size_t>(e));
    }
    listCuts();
    graph_.copy = CutCopy::cut(mesh_, graph_.cutEdges, graph_.cutEdgeBits);
    progress_.check();

    const BitArray rest = DiskCut::restCuts(
        {graph_.copy, adjacency_, graph_.mirrorHalf, field_.singularities.points, graph_.curveEdges}, progress_);
    progress_.check();
    if (!rest.any()) {
        if (symmetric_ == nullptr) {
            return;
        }
    } else {
        addDiskCuts(rest);
    }
    if (symmetric_ != nullptr) {
        mirrorCuts();
    }
    graph_.copy = CutCopy::cut(mesh_, graph_.cutEdges, graph_.cutEdgeBits);
}

BitArray Computation::treeEdges(std::vector<int>& planeChainEdges)
{
    // Curve edges are cheap, edges of the mirrored half a little dearer.
    std::vector<float> weights = edgeLengths(mesh_);
    for (std::size_t e = 0; e < edgeCount_; ++e) {
        if (graph_.curveEdgeKind[e] != curve_edge::none) {
            weights[e] = weights[e] / 3.0f;
        }
    }
    if (symmetric_ != nullptr) {
        forEachEdgeOfTriangles(adjacency_, graph_.mirrorHalf, [&](int e) {
            weights[static_cast<std::size_t>(e)] = weights[static_cast<std::size_t>(e)] * 1.1f;
        });
    }

    SteinerTree::Sources sources;
    sources.points = field_.singularities.points;
    sources.points.insert(sources.points.end(), graph_.curveCheckSources.begin(),
                          graph_.curveCheckSources.end());
    if (symmetric_ != nullptr) {
        // Every connected piece of the mirror plane is cut and is one source.
        const Csr chains = EdgeGroups::byConnectivity(adjacency_, symmetric_->planeEdges(adjacency_));
        for (int k = 0; k < chains.rows(); ++k) {
            if (!part_.mirror.classes.valid) {
                continue;
            }
            const auto chain = chains.row(k);
            planeChainEdges.insert(planeChainEdges.end(), chain.begin(), chain.end());
            std::vector<int> line;
            BitArray listed(static_cast<std::size_t>(mesh_.pointCount()));
            for (const int e : chain) {
                for (const int v : adjacency_.edgeVertices(e)) {
                    if (!listed.test(static_cast<std::size_t>(v))) {
                        listed.set(static_cast<std::size_t>(v));
                        line.push_back(v);
                    }
                }
            }
            sources.lines.push_back(std::move(line));
        }
    }
    // Border loops; with frozen borders only their kept points.
    const BoundaryLoops loops(mesh_);
    const auto* marks = prep::FreezeMarks::find(mesh_);
    for (int r = 0; r < loops.count(); ++r) {
        const auto loop = loops.points().row(r);
        std::vector<int> line;
        for (const int v : loop) {
            if (marks == nullptr || marks->values()[static_cast<std::size_t>(v)].value == 1) {
                line.push_back(v);
            }
        }
        sources.lines.push_back(std::move(line));
    }
    return SteinerTree::edges(adjacency_, weights, sources);
}

void Computation::addStars()
{
    // Around a singular point with a whole turn every edge is cut whose other
    // end is not on a cut yet.
    const Singularities& singular = field_.singularities;
    BitArray& bits = graph_.cutEdgeBits;
    for (std::size_t i = 0; i < singular.points.size(); ++i) {
        if ((singular.fullTurn[i] & 1) == 0) {
            continue;
        }
        const int v = singular.points[i];
        for (const int e : adjacency_.pointEdges(v)) {
            if (adjacency_.edgeCorners(e).size() <= 1) {
                continue;
            }
            const auto atOther = adjacency_.pointEdges(adjacency_.otherEnd(e, v));
            if (std::any_of(atOther.begin(), atOther.end(),
                            [&](int other) { return bits.test(static_cast<std::size_t>(other)); })) {
                continue;
            }
            bits.set(static_cast<std::size_t>(e));
        }
    }
}

void Computation::addDiskCuts(const BitArray& rest)
{
    // The cuts were found on the copy; they are kept as edges of the part mesh.
    std::vector<int> cuts;
    rest.forEachSet([&](std::size_t e) {
        for (const int original : graph_.copy.originalEdges(adjacency_, static_cast<int>(e))) {
            if (std::find(cuts.begin(), cuts.end(), original) == cuts.end()) {
                cuts.push_back(original);
            }
        }
    });
    if (!cuts.empty()) {
        graph_.diskCutBits = BitArray(edgeCount_);
    }
    for (const int e : cuts) {
        graph_.diskCutBits.set(static_cast<std::size_t>(e));
        graph_.cutEdgeBits.set(static_cast<std::size_t>(e));
    }
    graph_.diskCuts = std::move(cuts);
    listCuts();
    progress_.check();
}

void Computation::mirrorCuts()
{
    // The mirrored half takes the cuts of the computed half: cuts inside it
    // go, the mirror edges of the other cuts come.
    const std::vector<int>& edgeMirror = symmetric_->edgeMirror(adjacency_);
    const BitArray& half = graph_.mirrorHalf;
    auto insideHalf = [&](int e) {
        const auto faces = adjacency_.innerEdgeFaces(e);
        return faces && half.testBounded(static_cast<std::size_t>((*faces)[0])) &&
               half.testBounded(static_cast<std::size_t>((*faces)[1]));
    };
    std::vector<int> mirrored;
    BitArray& bits = graph_.cutEdgeBits;
    bits.forEachSet([&](std::size_t next) {
        const int e = static_cast<int>(next);
        if (adjacency_.edgeCorners(e).size() != 2) {
            return;
        }
        if (insideHalf(e)) {
            bits.set(next, false);
            return;
        }
        const int m = next < edgeMirror.size() ? edgeMirror[next] : -1;
        if (m >= 0 && insideHalf(m)) {
            mirrored.push_back(m);
        }
    });
    for (const int e : mirrored) {
        bits.set(static_cast<std::size_t>(e));
    }
    listCuts();
}

} // namespace

void CutGraphBuilder::compute(parts::Part& part, const Options& options, const run::Progress& progress)
{
    Computation(part, options, progress).run();
}

} // namespace areko::remesher::field
