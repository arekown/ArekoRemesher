#include "part_edge_length.h"

#include "prep/freeze_marks.h"
#include "prep/reference_length.h"

#include <cmath>

namespace areko::remesher::field {

void PartEdgeLength::assign(parts::Part& part, int quadCount, std::size_t partCount, float borderDensityInverse)
{
    const float area = part.mesh.triangleGeometry().totalArea();
    part.bordersFrozen = false;
    if (const auto* marks = prep::FreezeMarks::find(part.mesh)) {
        const mesh::BitArray* locked = part.lockedBorderEdges.empty() ? nullptr : &part.lockedBorderEdges;
        const float chord =
            prep::ReferenceLength::borderLength(part.mesh, &marks->values(), locked, borderDensityInverse);
        if (chord > 0.0f) {
            part.edgeLength = chord;
            part.bordersFrozen = true;
            return;
        }
    }
    if (!(area > 0.0f)) {
        return;
    }
    int n = part.quadBudget;
    if (n < 20) {
        if (quadCount > 400 || partCount >= 3) {
            n = 20;
        } else {
            const int least = part.mesh.adjacency().boundaryEdges().empty() ? 6 : 4;
            n = n < least ? least : n;
        }
    }
    const float given = part.targetEdgeLength;
    part.edgeLength = given > 0.0f ? given : std::sqrt(area / static_cast<float>(n));
}

} // namespace areko::remesher::field
