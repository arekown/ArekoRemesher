#pragma once

#include "param/constraint.h"
#include "param/system_context.h"

#include <vector>

namespace areko::remesher::param {

// Energy of the parameterisation: four rows per triangle over the variables
// of its three copy points. With the triangle laid flat (first corner at the
// origin, second on the first frame axis) and its gradients scaled by the
// target sizes h_u, h_v of the triangle, they ask that
//   - u grows by 1 per h_u along the first combed direction,
//   - u does not change along the second,
//   - v grows by 1 per h_v along the second,
//   - v does not change along the first.
// Twice the area of the flat triangle is replaced by the length of the cross
// product of its edges in space; the rows that ask for a growth of 1 are
// weighted by that length times a fixed factor, the others by the length.
// Triangles of the mirrored half get no rows when symmetry is switched on.
class EnergyRows {
public:
    struct Result {
        std::vector<TermList> rows;        // over variables, sorted by variable
        std::vector<double> targets;
        std::vector<double> weights;
        std::vector<int> triangleRowStart; // triangle -> its first row, one entry more than triangles
    };
    [[nodiscard]] static Result build(const SystemContext& context);
};

} // namespace areko::remesher::param
