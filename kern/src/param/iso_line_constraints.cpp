#include "iso_line_constraints.h"

#include "mesh/bit_array.h"
#include "param/variable_flags.h"

namespace areko::remesher::param {

using namespace areko::remesher::mesh;

namespace {

// The points of a run of edges of the cut copy: per edge its smaller end
// point, then its larger one, each point once at its first appearance.
[[nodiscard]] std::vector<int> pointsAlong(const Adjacency& copy, std::span<const int> edges)
{
    BitArray seen(static_cast<std::size_t>(copy.pointCount()));
    std::vector<int> points;
    for (const int e : edges) {
        for (const int w : copy.edgeVertices(e)) {
            if (!seen.test(static_cast<std::size_t>(w))) {
                seen.set(static_cast<std::size_t>(w));
                points.push_back(w);
            }
        }
    }
    return points;
}

} // namespace

IsoLineConstraints::Result IsoLineConstraints::build(const SystemContext& context, System& system)
{
    const CurveChains& curves = context.chains().curves;
    const Adjacency& copy = context.copy().adjacency();
    ConstraintSet& set = system.constraints;
    const std::vector<Constraint>& conditions = set.list();
    // Conditions from here on are the new ones.
    const int base = set.count();
    Result out;
    system.curveReference.assign(static_cast<std::size_t>(curves.count()), -1);

    for (int chain = 0; chain < curves.count(); ++chain) {
        const int coordinate = curves.coordinate[static_cast<std::size_t>(chain)];
        const auto edges = curves.edges.row(chain);
        if (coordinate < 0 || edges.empty()) {
            continue;
        }
        const std::vector<int> points = pointsAlong(copy, edges);
        int reference = -1;
        for (const int w : points) {
            if (set.ofVariable(coordinate + 2 * w) == -1) {
                reference = coordinate + 2 * w;
                break;
            }
        }
        system.curveReference[static_cast<std::size_t>(chain)] = reference;
        int referenceCondition = -1;
        if (reference == -1) {
            reference = coordinate + 2 * points.front();
            referenceCondition = set.ofVariable(reference);
        }
        const bool whole = curves.edgeKind[static_cast<std::size_t>(edges.front())] == field::curve_edge::crease;
        if (whole) {
            VariableFlags::markIfContinuous(system, reference, variable_flag::node);
        }
        for (const int w : points) {
            const int variable = coordinate + 2 * w;
            if (variable == reference) {
                continue;
            }
            if (whole) {
                VariableFlags::markIfContinuous(system, variable, variable_flag::node);
            }
            const int condition = set.ofVariable(variable);
            if (condition == -1) {
                set.add(variable, {Term{reference, 1.0}});
                continue;
            }
            // Expression of the variable minus expression of the reference.
            TermList row = condition < base ? conditions[static_cast<std::size_t>(condition)].terms
                                            : TermList{Term{variable, 1.0}};
            if (referenceCondition != -1) {
                for (const Term& term : conditions[static_cast<std::size_t>(referenceCondition)].terms) {
                    row.push_back({term.index, -term.value});
                }
            } else {
                row.push_back({reference, -1.0});
            }
            out.rows.push_back(std::move(row));
            if (out.preferred.empty() || out.preferred.back() != reference) {
                out.preferred.push_back(reference);
            }
        }
    }
    set.resolveAmong(base);
    return out;
}

} // namespace areko::remesher::param
