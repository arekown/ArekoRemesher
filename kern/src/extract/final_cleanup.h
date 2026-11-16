#pragma once

#include "raw_mesh.h"

namespace areko::remesher::extract {

// The last cleanup of the extraction: degenerate faces and corners go (the
// cleanup of the raw mesh, with every map kept: point maps such as the
// colours stay, faces, corners and half edges follow their origin and the
// hard edges the half edges), then every point without corner. The point
// data of the raw mesh follow the new numbering: per point lists keep the
// entries of the points that stay, neighbour numbers and the output point of
// every point of the triangle mesh are renumbered (-1 for a removed point).
// The edges along curves are dropped (their numbers are no longer valid).
class FinalCleanup {
public:
    static void run(RawMesh& raw);
};

} // namespace areko::remesher::extract
