#pragma once

#include "output_points.h"

#include "mesh/attribute_map.h"
#include "mesh/bit_array.h"
#include "mesh/boundary_loops.h"
#include "mesh/csr.h"
#include "mesh/mesh.h"

#include <span>
#include <vector>

namespace areko::remesher::extract {

// What frozen borders add.
struct FrozenBorderInput {
    const std::vector<mesh::FreezeMark>& marks;   // per point of the triangle mesh; 1 = kept border point
    // Point of the triangle mesh -> output point (-1 for none); every kept
    // border point has one, the gap rings read it without a check.
    const std::vector<int>& originalToOutput;
    const mesh::BitArray* lockedEdges = nullptr;  // locked border edges of the triangle mesh; null or empty: none
};

// Closes the raw mesh along its borders:
//  - points with several face fans among the border points that lie on no
//    border loop, and the end points of edges with more than two corners, get
//    one copy per further fan; the point data follow,
//  - with rebuilt polygons every border loop below 1000 points that reaches
//    no border of the triangle mesh becomes a polygon against the loop
//    direction, a four point loop closes instead by joining two opposite
//    points when that diagonal is below a ninth of the other and below
//    area / (500 quad budget),
//  - without frozen borders and with followed borders, the stretches of a
//    loop between two points on border edges become polygons,
//  - with frozen borders the rings between the loops that reach the border
//    and the kept border chains become polygons, cut at points on locked
//    border edges.
// New polygons are appended in this order; joined points take over the
// corners of the other point. The faces are set directly as polygons.
// One object closes once: run() adds to what it collected before.
class HoleClosing {
public:
    struct Options {
        bool rebuildPolygons = true;
        bool followBorders = true;
    };

    HoleClosing(mesh::Mesh& output, OutputPoints& points, const mesh::Mesh& original, int quadBudget,
                const FrozenBorderInput* frozen, Options options);

    void run();

private:
    void splitFans();
    void closeLoops(const mesh::BoundaryLoops& loops);
    [[nodiscard]] bool collapse(std::span<const int> quad, float threshold);
    void closeNotches(const mesh::BoundaryLoops& loops);
    void closeGapRings();
    [[nodiscard]] bool marked(int point) const;
    void addPolygon(std::span<const int> polygon);
    void apply();

    mesh::Mesh& output_;
    OutputPoints& points_;
    const mesh::Mesh& original_;
    int quadBudget_;
    const FrozenBorderInput* frozen_;
    Options options_;

    std::vector<int> polygonEnds_; // end of every new polygon in polygonCorners_
    std::vector<int> polygonCorners_;
    std::vector<int> joinedCorners_;
    std::vector<int> joinTargets_;
    mesh::Csr reachingLoops_; // loops that reach the border, kept for the gap rings
};

} // namespace areko::remesher::extract
