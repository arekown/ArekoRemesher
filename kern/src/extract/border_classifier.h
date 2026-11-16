#pragma once

#include "raw_mesh.h"

#include <span>
#include <vector>

namespace areko::remesher::parts {
struct Part;
} // namespace areko::remesher::parts

namespace areko::remesher::extract {

// Classifies the border of a raw mesh:
//  - along every border loop the points on border points of the triangle
//    mesh are anchors; between two anchors joined by a border edge of the
//    triangle mesh, every point not yet on a border edge that lies closer
//    than 0.0001 (absolute) to that edge is marked as on a border edge,
//  - with frozen borders, the edges inside the mesh at every point on a
//    border edge that does not come from a point of the triangle mesh (and
//    not from a locked border edge) are dissolved, with nearly straight
//    points.
class BorderClassifier {
public:
    static void run(RawMesh& raw, const parts::Part& part);

    // The points of one border loop in the order they are looked at: from
    // the first anchor round the loop and back to it, so that the stretch
    // before the first anchor is checked as well. A loop without anchor as
    // it is.
    [[nodiscard]] static std::vector<int> visitOrder(std::span<const int> loop, const RawMesh& raw);

private:
    [[nodiscard]] static bool anchor(const RawMesh& raw, int point);
};

} // namespace areko::remesher::extract
