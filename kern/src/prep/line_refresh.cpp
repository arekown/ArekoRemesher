#include "line_refresh.h"

#include "flat_triangle_check.h"
#include "mesh/measures.h"
#include "sliver_fixer.h"

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr int flipRuns = 5;
constexpr int flipMode = 3;

} // namespace

void LineRefresh::run(Mesh& mesh, PrepState& state)
{
    state.lines.extractFromSharpEdges(mesh);
    if (mesh.topology().mode() != FaceMode::Triangles) {
        return;
    }
    // First protected set: group borders with the mesh border (or the border
    // edges) and the edges of all lines.
    BitArray protectedEdges = state.lines.collectEdges(mesh, true, false, true);
    if (state.lines.state() == 1) {
        NetLines::addEdges(mesh, state.lines.lines(), false, protectedEdges);
    }
    for (int round = 0; round < flipRuns; ++round) {
        FlatTriangleBlock block;
        block.mode = flipMode;
        block.round = round;
        FlatTriangleCheck::run(mesh, protectedEdges, block);
        if (block.count > 0) {
            protectedEdges =
                state.lines.collectEdges(mesh, PrepState::useGroupBorders, state.hasEdgeAttribute, true);
        }
        if (!block.again) {
            break;
        }
    }
    BitArray protectedPoints;
    addEdgeEndPoints(mesh.adjacency(), protectedEdges, mesh.pointCount(), protectedPoints);
    SliverFixer::run(mesh, protectedPoints);
}

} // namespace areko::remesher::prep
