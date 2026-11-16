#include "anchor_constraints.h"

#include "mesh/boundary_loops.h"
#include "mesh/connected_parts.h"

#include <optional>
#include <span>
#include <vector>

namespace areko::remesher::param {

using namespace areko::remesher::mesh;

namespace {

// Border loops of the cut copy, piece by piece, from its neighbourhood and
// its pieces. They are only needed for a piece without singular point and
// are built when first asked for.
class CopyBorders {
public:
    CopyBorders(const field::CutCopy& copy, const ConnectedParts& pieces) : copy_(copy), pieces_(pieces) {}

    // The loops of piece `piece` in list order.
    [[nodiscard]] const std::vector<int>& loopsOf(int piece)
    {
        build();
        return loopsOfPiece_[static_cast<std::size_t>(piece)];
    }
    [[nodiscard]] std::span<const int> points(int loop) const { return loops_->points().row(loop); }

private:
    void build()
    {
        if (loops_) {
            return;
        }
        const Adjacency& adjacency = copy_.adjacency();
        loops_.emplace(adjacency, pieces_);
        loopsOfPiece_.assign(static_cast<std::size_t>(pieces_.partCount()), std::vector<int>{});
        for (int loop = 0; loop < loops_->count(); ++loop) {
            const int corner = adjacency.edgeCorners(loops_->edges().row(loop)[0])[0];
            const int piece = pieces_.partOfFace(adjacency.cornerFace(corner));
            loopsOfPiece_[static_cast<std::size_t>(piece)].push_back(loop);
        }
    }

    const field::CutCopy& copy_;
    const ConnectedParts& pieces_;
    std::optional<BoundaryLoops> loops_;
    std::vector<std::vector<int>> loopsOfPiece_;
};

// Whether curve edges of both combed directions meet at copy point `w`.
[[nodiscard]] bool curvesCross(const SystemContext& context, int w)
{
    const field::CutGraph& graph = context.graph();
    const Adjacency& adjacency = context.adjacency();
    const int origin = context.copy().pointOrigins()[static_cast<std::size_t>(w)];
    bool alongFirst = false;
    bool alongSecond = false;
    for (const int e : adjacency.pointEdges(origin)) {
        if (graph.curveEdgeKind[static_cast<std::size_t>(e)] == field::curve_edge::none) {
            continue;
        }
        const int corner = adjacency.edgeCorners(e)[0];
        if (graph.curveFollowsU.testBounded(static_cast<std::size_t>(corner))) {
            alongFirst = true;
        } else {
            alongSecond = true;
        }
    }
    return alongFirst && alongSecond;
}

// Anchor of a piece without singular point.
[[nodiscard]] int borderAnchor(const SystemContext& context, const ConnectedParts& pieces, CopyBorders& borders,
                               int piece)
{
    const std::vector<int>& loops = borders.loopsOf(piece);
    for (const int loop : loops) {
        for (const int w : borders.points(loop)) {
            if (!context.mirrored(w) && curvesCross(context, w)) {
                return w;
            }
        }
    }
    for (const int loop : loops) {
        for (const int w : borders.points(loop)) {
            if (!context.mirrored(w)) {
                return w;
            }
        }
    }
    return context.copyPoint(3 * pieces.partFaces().row(piece)[0]);
}

} // namespace

void AnchorConstraints::add(const SystemContext& context, System& system)
{
    const field::CutCopy& copy = context.copy();
    const ConnectedParts pieces(copy.topology(), copy.adjacency());
    const int pieceCount = pieces.partCount();
    std::vector<int> anchors(static_cast<std::size_t>(pieceCount), -1);
    int found = 0;
    for (const int singular : context.chains().singularities.points) {
        if (found == pieceCount) {
            break;
        }
        for (const int w : context.copies(singular)) {
            const int piece = pieces.partOfPoint(w);
            if (piece >= 0 && anchors[static_cast<std::size_t>(piece)] == -1) {
                anchors[static_cast<std::size_t>(piece)] = w;
                if (++found == pieceCount) {
                    break;
                }
            }
        }
    }
    CopyBorders borders(copy, pieces);
    for (int piece = 0; piece < pieceCount; ++piece) {
        int& anchor = anchors[static_cast<std::size_t>(piece)];
        if (anchor == -1) {
            anchor = borderAnchor(context, pieces, borders, piece);
        }
        system.constraints.add(2 * anchor, {});
        system.constraints.add(2 * anchor + 1, {});
    }
    system.constraints.closeAnchors();
    system.sections.anchors = system.constraints.count();
}

} // namespace areko::remesher::param
