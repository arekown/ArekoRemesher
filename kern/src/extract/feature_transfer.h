#pragma once

namespace areko::remesher::parts {
struct Part;
} // namespace areko::remesher::parts

namespace areko::remesher::extract {

// Brings the features of the triangle mesh of a part onto its quad mesh:
//  - without crease curves every face takes the face group of the first
//    triangle (when there are groups), nothing else,
//  - with crease curves the output points on them are marked (they give
//    their weight to their neighbours), the marked points are sorted into one
//    polyline per crease chain by their place along it, the chains are cut
//    into the quad mesh (or, without cutting, the edges between consecutive
//    polyline points become curve edges), then the sharp edges and the face
//    groups follow (HardEdgeTransfer, GroupTransfer).
class FeatureTransfer {
public:
    static void run(parts::Part& part, bool cutAlongCurves);
};

} // namespace areko::remesher::extract
