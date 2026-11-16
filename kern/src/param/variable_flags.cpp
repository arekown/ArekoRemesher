#include "variable_flags.h"

#include "mesh/boundary_loops.h"

#include <algorithm>

namespace areko::remesher::param {

using namespace areko::remesher::mesh;

namespace {

void mark(System& system, int variable, std::uint8_t flag)
{
    system.flags[static_cast<std::size_t>(variable)] = flag;
    system.periods[static_cast<std::size_t>(variable)] = 1;
}

[[nodiscard]] std::uint8_t flagOf(const System& system, int variable)
{
    return system.flags[static_cast<std::size_t>(variable)];
}

// Marks both coordinates of copy point `w` as kept free if both are marked.
bool keepFree(System& system, int w)
{
    if (flagOf(system, 2 * w) == variable_flag::continuous || flagOf(system, 2 * w + 1) == variable_flag::continuous) {
        return false;
    }
    mark(system, 2 * w, variable_flag::keptFree);
    mark(system, 2 * w + 1, variable_flag::keptFree);
    return true;
}

// Border points of a frozen border that stay where they are: every one but
// those marked as lying between kept points.
[[nodiscard]] bool roundedOnBorder(FreezeMark markOfPoint) noexcept
{
    return markOfPoint.value != 2 && markOfPoint.value != 3;
}

} // namespace

void VariableFlags::markIfContinuous(System& system, int variable, std::uint8_t flag)
{
    if (flagOf(system, variable) == variable_flag::continuous) {
        mark(system, variable, flag);
    }
}

void VariableFlags::assign(const SystemContext& context, System& system)
{
    const Chains& chains = context.chains();
    const Adjacency& adjacency = context.adjacency();
    system.copyPointCount = context.copy().pointCount();
    system.variableCount = 2 * system.copyPointCount + 2 * static_cast<int>(chains.genusChains.size());
    system.flags.assign(static_cast<std::size_t>(system.variableCount), variable_flag::continuous);
    system.periods.assign(static_cast<std::size_t>(system.variableCount), 0);
    system.periodChanges = 0;
    system.constraints = ConstraintSet(system.variableCount);

    for (std::size_t i = 0; i < chains.nodes.size(); ++i) {
        const int node = chains.nodes[i];
        const auto edges = adjacency.pointEdges(node);
        const bool onBorder =
            std::any_of(edges.begin(), edges.end(), [&](int e) { return adjacency.isBoundaryEdge(e); });
        if (!chains.singularNodes.testBounded(i) && !onBorder) {
            continue;
        }
        for (const int c : context.copies(node)) {
            mark(system, 2 * c, variable_flag::node);
            mark(system, 2 * c + 1, variable_flag::node);
        }
    }
    for (int variable = 2 * system.copyPointCount; variable < system.variableCount; ++variable) {
        mark(system, variable, variable_flag::offset);
    }
    for (const int source : context.graph().curveCheckSources) {
        for (const int c : context.copies(source)) {
            markIfContinuous(system, 2 * c, variable_flag::checkSource);
            markIfContinuous(system, 2 * c + 1, variable_flag::checkSource);
        }
    }
    if (const std::vector<FreezeMark>* marks = context.part().frozenBorderMarks()) {
        const BoundaryLoops loops(context.mesh());
        for (const int p : loops.points().entries) {
            if (!roundedOnBorder((*marks)[static_cast<std::size_t>(p)])) {
                continue;
            }
            for (const int c : context.copies(p)) {
                if (!context.mirrored(c)) {
                    markIfContinuous(system, 2 * c, variable_flag::frozenBorder);
                    markIfContinuous(system, 2 * c + 1, variable_flag::frozenBorder);
                }
            }
        }
    }
}

std::vector<int> VariableFlags::markChainEnds(const SystemContext& context, System& system)
{
    const Chains& chains = context.chains();
    const field::CutGraph& graph = context.graph();
    const Adjacency& adjacency = context.adjacency();
    std::vector<int> preferred;
    for (const int node : chains.nodes) {
        for (const int c : context.copies(node)) {
            for (const int variable : {2 * c, 2 * c + 1}) {
                if (flagOf(system, variable) == variable_flag::continuous) {
                    preferred.push_back(variable);
                }
            }
        }
    }
    if (!graph.diskCuts.empty()) {
        for (int chain = 0; chain < chains.count(); ++chain) {
            const auto edges = chains.edges.row(chain);
            if (edges.empty() || chains.offsetVariable[static_cast<std::size_t>(chain)] != 0) {
                continue;
            }
            for (const int e : {edges.front(), edges.back()}) {
                if (!graph.diskCutBits.testBounded(static_cast<std::size_t>(e))) {
                    continue;
                }
                const auto corners = adjacency.edgeCorners(e);
                if (corners.size() != 2) {
                    continue;
                }
                // The copy points of the edge in its first and in its second triangle.
                bool first = keepFree(system, context.copyPoint(corners[0]));
                first = keepFree(system, context.copyPoint(adjacency.next(corners[0]))) || first;
                bool second = keepFree(system, context.copyPoint(corners[1]));
                second = keepFree(system, context.copyPoint(adjacency.next(corners[1]))) || second;
                if (first && second) {
                    break;
                }
            }
        }
    }
    for (const int node : chains.nodes) {
        const std::vector<int> copies = context.copies(node);
        for (std::size_t i = 1; i < copies.size(); ++i) {
            const int c = copies[i];
            if (flagOf(system, 2 * c) == variable_flag::keptFree) {
                continue;
            }
            for (const int variable : {2 * c, 2 * c + 1}) {
                if (flagOf(system, variable) != variable_flag::continuous) {
                    preferred.push_back(variable);
                    mark(system, variable, variable_flag::nodeCopy);
                }
            }
        }
    }
    return preferred;
}

} // namespace areko::remesher::param
