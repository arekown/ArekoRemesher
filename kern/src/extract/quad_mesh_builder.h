#pragma once

#include "raw_mesh.h"
#include "trace.h"

namespace areko::remesher::parts {
struct Part;
} // namespace areko::remesher::parts

namespace areko::remesher::extract {

// Builds the raw mesh of a part from its trace: one output point per
// numbered node (in number order), the faces closed by the trace as polygons
// with the merged nodes replaced by their targets, and per output point its
// source on the triangle mesh, its triangle, its neighbours along the four
// directions and whether it lies on or ran into the border.
//  - a node closer to a corner of its triangle than five millionths of the
//    mesh diagonal is moved onto that corner and comes from that point,
//  - with frozen borders, every kept border point without output point gets
//    one of its own after the nodes (unused by the faces).
class QuadMeshBuilder {
public:
    [[nodiscard]] static RawMesh build(const parts::Part& part, const Trace& trace);

    // Moves `node` onto the nearest corner of its triangle within the squared
    // distance `toleranceSq` (corners in order, each measured from the place
    // the node had; a later corner wins only when nearer): the node becomes
    // a node on that point. Seed nodes stay.
    static void snapToCorner(const mesh::Mesh& triangles, double toleranceSq, TraceNode& node);
};

} // namespace areko::remesher::extract
