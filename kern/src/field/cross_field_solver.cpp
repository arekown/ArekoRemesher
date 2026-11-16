#include "cross_field_solver.h"

#include "cross_field_rounding.h"
#include "cross_field_system.h"
#include "cross_frames.h"
#include "dual_spanning_tree.h"
#include "singularity_finder.h"

#include "solver/constrained_driver.h"

#include <optional>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

namespace {

// The working state of one cross field computation.
class Computation {
public:
    Computation(parts::Part& part, Curvature& curvature, const CrossFieldSolver::Options& options,
                const run::Progress& progress)
        : part_(part), mesh_(part.mesh), adjacency_(part.mesh.adjacency()), curvature_(curvature), options_(options),
          progress_(progress), faceCount_(static_cast<std::size_t>(part.mesh.faceCount())),
          edgeCount_(static_cast<std::size_t>(part.mesh.adjacency().edgeCount()))
    {
    }

    void run();

private:
    void findMirrorPairs();
    void selectAnchors();
    [[nodiscard]] bool solve();
    void mirrorAngles(BitArray& marked);
    void roundJumps();
    void matchMirrorJumps(BitArray& marked);
    // Where an edge lies against the mirror plane: 0 or 1 its side, 2 on or
    // across the plane, -1 undecided.
    [[nodiscard]] int edgeSide(int e) const;

    parts::Part& part_;
    const Mesh& mesh_;
    const Adjacency& adjacency_;
    Curvature& curvature_;
    CrossFieldSolver::Options options_;
    const run::Progress& progress_;
    std::size_t faceCount_;
    std::size_t edgeCount_;
    // Mirror pairs of a symmetric part, and the pairs that mirror every
    // triangle (those of a symmetric part, or what is left of older pairs).
    prep::SymmetryCache* symmetric_ = nullptr;
    prep::SymmetryCache* paired_ = nullptr;
    std::optional<prep::SymmetryCache> leftover_;
    BitArray anchors_;
    BitArray freeEdges_;
    std::vector<float> connection_;
    CrossField field_;
};

void Computation::run()
{
    part_.crossField.reset();
    findMirrorPairs();
    progress_.check();
    selectAnchors();
    progress_.check();
    // The mirrored half of the field is the half the curvature did not compute.
    freeEdges_ = DualSpanningTree::freeEdges(mesh_, anchors_, curvature_.mirrorFaces);
    progress_.check();
    CrossFrames::compute(mesh_, symmetric_, field_);
    progress_.check();
    connection_ = CrossFrames::connection(mesh_, field_.frameE1);
    progress_.check();
    if (!solve()) {
        return;
    }
    progress_.check();
    field_.singularities = SingularityFinder::find(mesh_, field_, connection_, symmetric_, part_.mirror.classes.valid);
    progress_.check();
    part_.crossField = std::move(field_);
}

void Computation::findMirrorPairs()
{
    symmetric_ = part_.mirror.evaluate(mesh_, options_.symmetry);
    if (!options_.symmetry) {
        return;
    }
    if (symmetric_ != nullptr) {
        paired_ = symmetric_->allFacesMirrored(mesh_) ? symmetric_ : nullptr;
    } else if (const prep::SymmetryCache* stored = part_.mirror.stored()) {
        leftover_.emplace(*stored);
        paired_ = leftover_->allFacesMirrored(mesh_) ? &*leftover_ : nullptr;
    }
}

void Computation::selectAnchors()
{
    anchors_ = BitArray(faceCount_);
    std::vector<float>& weight = curvature_.weight;
    bool any = false;
    for (std::size_t t = 0; t < faceCount_; ++t) {
        if (weight[t] >= CrossFieldSystem::anchorWeight) {
            anchors_.set(t);
            any = true;
        }
    }
    if (!any && faceCount_ > 0) {
        anchors_.set(0);
        weight[0] = CrossFieldSystem::anchorWeight;
    }
    if (curvature_.mirrorFaceCount == 0 || symmetric_ == nullptr) {
        return;
    }
    // The half that is not computed takes the curvature values of its mirror triangles.
    const std::vector<int>& faceMirror = symmetric_->faceMirror(mesh_);
    for (std::size_t t = 0; t < faceCount_ && t < faceMirror.size(); ++t) {
        const int mirrored = faceMirror[t];
        if (mirrored == -1 || !curvature_.mirrorFaces.testBounded(static_cast<std::size_t>(mirrored))) {
            continue;
        }
        const std::size_t m = static_cast<std::size_t>(mirrored);
        curvature_.lambda1[m] = curvature_.lambda1[t];
        curvature_.lambda2[m] = curvature_.lambda2[t];
        curvature_.flags[m] |= curvature_flags::reliable;
        weight[m] = weight[t];
    }
}

bool Computation::solve()
{
    field_.targetAngle = CrossFieldSystem::targetAngles(mesh_, curvature_, field_.frameE1,
                                                        symmetric_ != nullptr ? &symmetric_->faceMirror(mesh_) : nullptr);
    const CrossFieldSystem::Input input{mesh_, curvature_, anchors_, freeEdges_, connection_, field_.targetAngle};
    field_.periodJump = CrossFieldSystem::anchorJumps(input);
    field_.angle.assign(faceCount_, 0.0f);
    progress_.check();

    CrossFieldSystem system(input, field_.periodJump);
    // The driver only solves; every integer is rounded by the rounding rounds.
    solver::DriverOptions driverOptions;
    driverOptions.keepGoing = true;
    driverOptions.localUpdate = true;
    driverOptions.local.tolerance = wide(0.1f);
    driverOptions.local.maxIterations = 50;
    driverOptions.fullSolveEvery = 250;
    CrossFieldRounding rounding(system, progress_, driverOptions);
    std::vector<double> x;
    if (!rounding.run(x)) {
        return false;
    }
    const std::vector<int>& faces = system.unknownFaces();
    for (std::size_t u = 0; u < faces.size(); ++u) {
        field_.angle[static_cast<std::size_t>(faces[u])] = static_cast<float>(x[u]);
    }
    BitArray marked(edgeCount_);
    mirrorAngles(marked);
    roundJumps();
    matchMirrorJumps(marked);
    return true;
}

void Computation::mirrorAngles(BitArray& marked)
{
    if (paired_ == nullptr) {
        return;
    }
    // The mirrored half takes pi minus the angle of its mirror triangle; its
    // edges on side 0 are matched to their mirror edges afterwards.
    const std::vector<int>& faceMirror = paired_->faceMirror(mesh_);
    const BitArray& half = curvature_.mirrorFaces;
    for (std::ptrdiff_t next = half.findNext(0); next != -1; next = half.findNext(static_cast<std::size_t>(next) + 1)) {
        const std::size_t t = static_cast<std::size_t>(next);
        field_.angle[t] = static_cast<float>(piD - wide(field_.angle[static_cast<std::size_t>(faceMirror[t])]));
        for (int i = 0; i < 3; ++i) {
            const int e = adjacency_.cornerEdge(3 * static_cast<int>(t) + i);
            if (e >= 0 && edgeSide(e) == 0) {
                marked.set(static_cast<std::size_t>(e));
            }
        }
    }
}

void Computation::roundJumps()
{
    // Quarter turns between the crosses as they came out, across every inner
    // edge that has a computed triangle.
    const BitArray& half = curvature_.mirrorFaces;
    for (std::size_t e = 0; e < edgeCount_; ++e) {
        const auto faces = adjacency_.innerEdgeFaces(static_cast<int>(e));
        if (!faces) {
            continue;
        }
        const std::size_t f = static_cast<std::size_t>((*faces)[0]);
        const std::size_t g = static_cast<std::size_t>((*faces)[1]);
        if (half.testBounded(f) && half.testBounded(g)) {
            continue;
        }
        field_.periodJump[e] = nearestQuarterTurns(wide(field_.angle[f]), wide(field_.angle[g]), connection_[e]);
    }
}

void Computation::matchMirrorJumps(BitArray& marked)
{
    if (paired_ == nullptr) {
        return;
    }
    const std::vector<int>& faceMirror = paired_->faceMirror(mesh_);
    const std::vector<int>& edgeMirror = paired_->edgeMirror(adjacency_);
    for (std::size_t e = 0; e < edgeCount_; ++e) {
        if (edgeSide(static_cast<int>(e)) != 0) {
            continue;
        }
        const auto faces = adjacency_.innerEdgeFaces(static_cast<int>(e));
        if (faces && anchors_.testBounded(static_cast<std::size_t>((*faces)[0])) &&
            anchors_.testBounded(static_cast<std::size_t>((*faces)[1]))) {
            marked.set(e);
        }
    }
    // A marked edge takes the quarter turn that fits its mirror edge modulo 4
    // and lies nearest to the turn of the field across it.
    for (std::ptrdiff_t next = marked.findNext(0); next != -1; next = marked.findNext(static_cast<std::size_t>(next) + 1)) {
        const std::size_t e = static_cast<std::size_t>(next);
        const int mirrored = e < edgeMirror.size() ? edgeMirror[e] : -1;
        const auto faces = adjacency_.innerEdgeFaces(static_cast<int>(e));
        const auto mirrorFaces = mirrored >= 0 ? adjacency_.innerEdgeFaces(mirrored) : std::nullopt;
        if (!faces || !mirrorFaces) {
            continue;
        }
        const std::size_t f = static_cast<std::size_t>((*faces)[0]);
        const std::size_t g = static_cast<std::size_t>((*faces)[1]);
        int fit = field_.periodJump[static_cast<std::size_t>(mirrored)];
        if (f < faceMirror.size() && (*mirrorFaces)[0] == faceMirror[f]) {
            fit = -fit;
        }
        const float turns =
            (((connection_[e] - field_.angle[f]) + field_.angle[g]) / quarterTurnF - static_cast<float>(fit)) * 0.25f;
        field_.periodJump[e] =
            static_cast<int>(static_cast<unsigned>(util::roundHalfAway(turns)) * 4u + static_cast<unsigned>(fit));
    }
}

int Computation::edgeSide(int e) const
{
    const std::vector<int>& pairs = paired_->pointMirror();
    const BitArray& side = paired_->pointSide();
    const auto ends = adjacency_.edgeVertices(e);
    const std::size_t u = static_cast<std::size_t>(ends[0]);
    const std::size_t v = static_cast<std::size_t>(ends[1]);
    if (u >= pairs.size() || v >= pairs.size()) {
        return -1;
    }
    const int mirrorU = pairs[u];
    const int mirrorV = pairs[v];
    const bool sideU = side.testBounded(u);
    const bool sideV = side.testBounded(v);
    if (mirrorU == ends[0]) {
        return mirrorV == ends[1] ? 2 : (sideV ? 1 : 0);
    }
    if (mirrorV == ends[1]) {
        return sideU ? 1 : 0;
    }
    if (mirrorU == -1 || mirrorV == -1) {
        return -1;
    }
    if (sideU == sideV) {
        return sideU ? 1 : 0;
    }
    return ends[1] == mirrorU ? 2 : -1;
}

} // namespace

void CrossFieldSolver::compute(parts::Part& part, const Options& options, const run::Progress& progress)
{
    if (!part.curvature) {
        part.crossField.reset();
        return;
    }
    Computation(part, *part.curvature, options, progress).run();
}

} // namespace areko::remesher::field
