#pragma once

#include "symmetry_cache.h"

#include "mesh/bit_array.h"
#include "mesh/mesh.h"

namespace areko::remesher::prep {

// Mirror evaluation of a part mesh against its mirror pairs; component 0 is
// the connected part of face 0.
class PartSymmetry {
public:
    // Valid when every face of component 0 has a mirror face and all mirror
    // faces lie in one component (`mirrorGroup`).
    struct Classes {
        bool valid = false;
        int mirrorGroup = -1;
    };

    [[nodiscard]] static Classes classify(const mesh::Mesh& mesh, SymmetryCache& cache);

    // Faces of component 0 wholly on side 0, when component 0 mirrors onto
    // itself and exactly half of its faces lie wholly on side 0 and the other
    // half wholly on side 1; otherwise no face. One bit per face.
    [[nodiscard]] static mesh::BitArray half(const mesh::Mesh& mesh, const SymmetryCache& cache,
                                             const Classes& classes);
};

} // namespace areko::remesher::prep
