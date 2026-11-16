#pragma once

#include "mesh/topology.h"
#include "param/system.h"
#include "param/system_context.h"
#include "prep/symmetry_cache.h"

#include <array>
#include <span>
#include <vector>

namespace areko::remesher::param {

// Turns a solution of a system (one value per unknown) into grid coordinates
// on the cut copy.
class UvOutput {
public:
    using Coordinates = std::vector<std::array<float, 2>>;

    // The value of every variable, in float: an unknown takes its value
    // rounded to float; an eliminated variable is the sum its condition
    // names, formed in double over the float values of the other variables
    // and rounded, the conditions taken from the last to the first. Every
    // other variable is 0.
    [[nodiscard]] static std::vector<float> variables(const System& system, std::span<const double> x);

    // (u, v) of every point of the cut copy for the solution `x`.
    [[nodiscard]] static Coordinates coordinates(const System& system, std::span<const double> x);

    // Gives the corners of the mirrored half the coordinates of their mirror
    // corners with u negated, triangle by triangle in ascending order. A
    // corner without mirror corner keeps its coordinates.
    static void mirror(const SystemContext& context, prep::SymmetryCache& pairs, Coordinates& uv);

    // The corners of the mirror triangle of triangle `t` that lie at the
    // mirror points of its three corners, in the order of its corners. The
    // mirror triangle runs the other way round; only the turn that its first
    // corner picks is tried. All -1 when the triangle has no mirror triangle
    // or its corners do not fit.
    [[nodiscard]] static std::array<int, 3> mirroredCorners(const mesh::Topology& topology,
                                                            std::span<const int> pointMirror,
                                                            std::span<const int> faceMirror, int t);
};

} // namespace areko::remesher::param
