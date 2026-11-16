#include "disk_cut.h"

#include "cut_graph.h"

#include "mesh/edge_chains.h"
#include "util/index_queue.h"
#include "util/reach_queue.h"

#include <vector>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

namespace {

// The working state of one search for the remaining cuts.
class Search {
public:
    Search(const DiskCut::Input& input, const run::Progress& progress)
        : input_(input), progress_(progress), copy_(input.copy.adjacency()),
          faceCount_(static_cast<std::size_t>(input.copy.topology().faceCount())),
          edgeCount_(static_cast<std::size_t>(input.copy.adjacency().edgeCount())), forest_(faceCount_)
    {
    }

    [[nodiscard]] BitArray run();

private:
    void markFreeEdges();
    void buildCurveChains();
    void growForest();
    [[nodiscard]] bool openLongestChain();
    void selectCandidates();
    void pruneDanglingEnds();
    void straighten();
    // Triangle across the edge of corner c, or -1.
    [[nodiscard]] int across(int c) const { return copy_.otherFace(copy_.cornerEdge(c), c / 3); }
    [[nodiscard]] bool isCandidate(int e) const { return e >= 0 && candidates_.test(static_cast<std::size_t>(e)); }

    const DiskCut::Input& input_;
    const run::Progress& progress_;
    const Adjacency& copy_;
    std::size_t faceCount_;
    std::size_t edgeCount_;
    BitArray free_;                 // edges that cost nothing
    BitArray blocked_;              // curve edges
    std::vector<int> chainOfEdge_;  // -1: none
    std::vector<int> chainLength_;  // edges per chain
    std::vector<int> parent_;       // -1: not reached; a root is its own parent
    util::ReachQueue forest_;       // search over the triangles; reach in crossed edges that cost
    std::vector<int> waiting_;      // corners whose curve edge has an unreached triangle behind it
    BitArray isWaiting_;
    BitArray candidates_;
    BitArray onBorder_;             // points of the copy with a border edge
    std::vector<int> degree_;       // candidate edges per point of the copy
};

BitArray Search::run()
{
    markFreeEdges();
    buildCurveChains();
    growForest();
    progress_.check();
    selectCandidates();
    pruneDanglingEnds();
    straighten();
    return std::move(candidates_);
}

void Search::markFreeEdges()
{
    free_ = BitArray(edgeCount_);
    BitArray singular(static_cast<std::size_t>(input_.original.pointCount()));
    for (const int v : input_.singularPoints) {
        singular.set(static_cast<std::size_t>(v));
    }
    const std::vector<int>& origin = input_.copy.pointOrigins();
    onBorder_ = BitArray(static_cast<std::size_t>(copy_.pointCount()));
    for (const int p : copy_.boundaryPoints()) {
        onBorder_.set(static_cast<std::size_t>(p));
        if (!singular.test(static_cast<std::size_t>(origin[static_cast<std::size_t>(p)]))) {
            for (const int e : copy_.pointEdges(p)) {
                free_.set(static_cast<std::size_t>(e));
            }
        }
    }
    forEachEdgeOfTriangles(copy_, input_.mirrorHalf, [&](int e) { free_.set(static_cast<std::size_t>(e)); });
}

void Search::buildCurveChains()
{
    blocked_ = BitArray(edgeCount_);
    for (const int edge : input_.curveEdges) {
        for (const int e : input_.copy.copyEdges(input_.original, edge)) {
            if (e >= 0) {
                blocked_.set(static_cast<std::size_t>(e));
            }
        }
    }
    Csr chainEdges;
    (void)EdgeChains::build(copy_, blocked_, &chainEdges);
    chainOfEdge_.assign(edgeCount_, -1);
    chainLength_.assign(static_cast<std::size_t>(chainEdges.rows()), 0);
    for (int k = 0; k < chainEdges.rows(); ++k) {
        chainLength_[static_cast<std::size_t>(k)] = chainEdges.rowSize(k);
        for (const int e : chainEdges.row(k)) {
            chainOfEdge_[static_cast<std::size_t>(e)] = k;
        }
    }
}

void Search::growForest()
{
    parent_.assign(faceCount_, -1);
    isWaiting_ = BitArray(3 * faceCount_);
    std::size_t root = 0;
    while (true) {
        while (root < faceCount_ && parent_[root] != -1) {
            ++root;
        }
        if (root == faceCount_) {
            return;
        }
        parent_[root] = static_cast<int>(root);
        forest_.start(static_cast<int>(root));
        do {
            forest_.run([&](int t) {
                for (int c = 3 * t; c < 3 * t + 3; ++c) {
                    const int e = copy_.cornerEdge(c);
                    if (e < 0 || copy_.edgeCorners(e).size() != 2) {
                        continue;
                    }
                    const int n = copy_.otherFace(e, t);
                    if (n < 0) {
                        continue;
                    }
                    const std::size_t ni = static_cast<std::size_t>(n);
                    if (free_.test(static_cast<std::size_t>(e))) {
                        // A free edge leads on at once.
                        if (forest_.offer(n, forest_.reach(t), true)) {
                            parent_[ni] = t;
                        }
                    } else if (blocked_.test(static_cast<std::size_t>(e))) {
                        if (parent_[ni] < 0 && !isWaiting_.test(static_cast<std::size_t>(c))) {
                            isWaiting_.set(static_cast<std::size_t>(c));
                            waiting_.push_back(c);
                        }
                    } else if (forest_.offer(n, forest_.reach(t) + 1.0f)) {
                        parent_[ni] = t;
                    }
                }
            });
        } while (openLongestChain());
    }
}

bool Search::openLongestChain()
{
    if (waiting_.empty()) {
        return false;
    }
    // Longest chain among the crossings that still lead to an unreached
    // triangle; the first one wins a tie.
    int best = -1;
    int longest = 0;
    for (const int c : waiting_) {
        const int chain = chainOfEdge_[static_cast<std::size_t>(copy_.cornerEdge(c))];
        if (parent_[static_cast<std::size_t>(across(c))] >= 0 || chain == -1 ||
            chainLength_[static_cast<std::size_t>(chain)] <= longest) {
            continue;
        }
        longest = chainLength_[static_cast<std::size_t>(chain)];
        best = chain;
    }
    std::vector<int> kept;
    if (best != -1) {
        // Exactly one crossing of that chain is opened: its first one with an
        // unreached triangle behind it. The crossings of the chain before
        // that one lapse; every other crossing, the later ones of the chain
        // included, waits on while the triangle behind it is unreached.
        int target = best;
        for (const int c : waiting_) {
            const std::size_t n = static_cast<std::size_t>(across(c));
            if (chainOfEdge_[static_cast<std::size_t>(copy_.cornerEdge(c))] == target) {
                if (parent_[n] < 0) {
                    // An unreached triangle takes any reach.
                    parent_[n] = c / 3;
                    (void)forest_.offer(static_cast<int>(n), forest_.reach(c / 3) + 1.0f);
                    target = -2;
                }
            } else if (parent_[n] == -1) {
                kept.push_back(c);
            }
        }
    }
    // Without a crossing left to open, the search goes on with the next root.
    for (const int c : waiting_) {
        isWaiting_.set(static_cast<std::size_t>(c), false);
    }
    for (const int c : kept) {
        isWaiting_.set(static_cast<std::size_t>(c));
    }
    waiting_ = std::move(kept);
    return best != -1;
}

void Search::selectCandidates()
{
    candidates_ = BitArray(edgeCount_, true);
    for (std::size_t f = 0; f < faceCount_; ++f) {
        // A root has no edge towards a parent.
        if (parent_[f] >= 0 && parent_[f] != static_cast<int>(f)) {
            const int e = copy_.commonEdge(static_cast<int>(f), parent_[f]);
            if (e != -1) {
                candidates_.set(static_cast<std::size_t>(e), false);
            }
        }
    }
    for (const int e : copy_.boundaryEdges()) {
        candidates_.set(static_cast<std::size_t>(e), false);
    }
    forEachEdgeOfTriangles(copy_, input_.mirrorHalf,
                           [&](int e) { candidates_.set(static_cast<std::size_t>(e), false); });
}

void Search::pruneDanglingEnds()
{
    const std::size_t pointCount = static_cast<std::size_t>(copy_.pointCount());
    degree_.assign(pointCount, 0);
    candidates_.forEachSet([&](std::size_t e) {
        for (const int p : copy_.edgeVertices(static_cast<int>(e))) {
            ++degree_[static_cast<std::size_t>(p)];
        }
    });
    util::IndexQueue ends(pointCount);
    for (std::size_t p = 0; p < pointCount; ++p) {
        if (degree_[p] == 1 && !onBorder_.test(p)) {
            ends.pushBack(static_cast<int>(p));
        }
    }
    while (!ends.empty()) {
        const int p = ends.pop();
        if (degree_[static_cast<std::size_t>(p)] == 0) {
            continue;
        }
        for (const int e : copy_.pointEdges(p)) {
            if (!candidates_.test(static_cast<std::size_t>(e))) {
                continue;
            }
            candidates_.set(static_cast<std::size_t>(e), false);
            degree_[static_cast<std::size_t>(p)] = 0;
            const std::size_t q = static_cast<std::size_t>(copy_.otherEnd(e, p));
            if (--degree_[q] == 1 && !onBorder_.test(q)) {
                ends.pushBack(static_cast<int>(q));
            }
            break;
        }
    }
}

void Search::straighten()
{
    // Two candidates of one triangle that meet in a point with no other
    // candidate give way to the third side of the triangle.
    auto replace = [&](int e, int with, int other, int x) {
        candidates_.set(static_cast<std::size_t>(other), false);
        candidates_.set(static_cast<std::size_t>(e), false);
        candidates_.set(static_cast<std::size_t>(with));
        degree_[static_cast<std::size_t>(x)] = 0;
    };
    auto lonely = [&](int x) {
        return degree_[static_cast<std::size_t>(x)] == 2 && !onBorder_.test(static_cast<std::size_t>(x));
    };
    for (int round = 0; round < 3; ++round) {
        bool changed = false;
        candidates_.forEachSet([&](std::size_t next) {
            const int e = static_cast<int>(next);
            for (const int c : copy_.edgeCorners(e)) {
                const int a = copy_.cornerEdge(copy_.prev(c));
                const int b = copy_.cornerEdge(copy_.next(c));
                if (a == -1 || b == -1) {
                    continue;
                }
                if (!isCandidate(a) && isCandidate(b) && !copy_.isBoundaryEdge(a) && lonely(copy_.vertex(copy_.next(c)))) {
                    replace(e, a, b, copy_.vertex(copy_.next(c)));
                    changed = true;
                    break;
                }
                if (isCandidate(a) && !isCandidate(b) && !copy_.isBoundaryEdge(b) && lonely(copy_.vertex(c))) {
                    replace(e, b, a, copy_.vertex(c));
                    changed = true;
                    break;
                }
            }
        });
        if (!changed) {
            break;
        }
    }
}

} // namespace

BitArray DiskCut::restCuts(const Input& input, const run::Progress& progress)
{
    return Search(input, progress).run();
}

} // namespace areko::remesher::field
