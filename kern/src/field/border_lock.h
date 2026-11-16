#pragma once

#include "mesh/vec3.h"
#include "parts/part.h"
#include "settings/settings.h"

namespace areko::remesher::field {

// Locks the border edges of a part that lie on the plane of a kept half: the
// cut left them as a border, but they are no border of the shape, so their
// length is not to shape the quads.
class BorderLock {
public:
    struct Options {
        bool halfMode = false;                // only one half of the shape is computed and kept
        bool axis[3] = {false, false, false}; // axes whose plane locks
        mesh::Vec3f planeMax;                 // largest coordinates of the prepared mesh
    };

    // The options for the symmetry switches as they are after the cut
    // (`planeMax` is left to the caller). The axes are the mirror axes. Only
    // a half that was handed in as such without naming an axis (keeping the
    // half is set, the symmetric computation is off) counts along all three.
    // A half mode without axis that comes from a guiding curve cutting into
    // a mesh whose symmetry is to be found has no axis: nothing is locked.
    [[nodiscard]] static Options optionsFor(const settings::SymmetrySwitches& symmetry);

    // Writes part.lockedBorderEdges: one bit per edge, set for a border edge
    // whose two end points have the same coordinate along one of the axes, and
    // that coordinate is the largest one of the prepared mesh (exact
    // comparison). Empty when no edge is locked.
    static void assign(parts::Part& part, const Options& options);
};

} // namespace areko::remesher::field
