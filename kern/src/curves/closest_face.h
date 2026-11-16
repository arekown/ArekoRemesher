#pragma once

#include "mesh/mesh.h"

#include <array>
#include <vector>

namespace areko::remesher::curves {

// Nearest triangle of a triangle mesh to a point within a radius. The point
// is dropped onto the plane of each triangle (its normal); inside the
// triangle (barycentric, float) that is the nearest point, otherwise the
// nearest point on the three sides (second side on ties with the first,
// third on ties with the better one). A triangle replaces the best one only
// when it is strictly closer.
//
// The triangles are visited in the order of a box tree, so among equally
// close triangles the one met first in that order wins:
//  - every node holds a box; the root the box of the mesh, every other node
//    the union of the boxes of its triangles,
//  - a node with at most 5 triangles is a leaf; otherwise it is split at the
//    middle of its box along its longest side (y before x, z before both on
//    equal sides): triangles whose box starts at or above the middle go
//    right, those ending at or below it go left (both in their order, the
//    right ones from the end backwards), and the ones across the middle
//    follow to the side they reach less far into,
//  - a split that leaves fewer than 2 triangles on one side makes a node of
//    at most 10 triangles a leaf, a larger one is cut at the middle index,
//  - the search walks the tree depth first, the right child before the
//    left, skips every node whose box the point is not within the best
//    distance of, and tries the triangles of a leaf in their order.
class ClosestFaceSearch {
public:
    explicit ClosestFaceSearch(const mesh::Mesh& mesh);

    struct Hit {
        int face = -1;
        mesh::Vec3f point;
        float distance = 0.0f;
    };
    // Face -1 when no triangle lies closer than `radius`.
    [[nodiscard]] Hit nearest(mesh::Vec3f p, float radius) const;

private:
    // Box as min x, y, z, max x, y, z.
    using Box = std::array<float, 6>;
    // A leaf holds the triangles `order_[begin, end)`, an inner node its two children.
    struct Node {
        Box box;
        bool leaf = false;
        int begin = 0;
        int end = 0;
        int left = -1;
        int right = -1;
    };

    // Nearest point of triangle f to p and its distance.
    [[nodiscard]] Hit onTriangle(int f, mesh::Vec3f p) const;
    int build(const Box& box, int begin, int end);
    [[nodiscard]] Box unionOf(int begin, int end) const;

    const mesh::Mesh& mesh_;
    std::vector<Box> faceBoxes_;
    std::vector<int> order_;
    std::vector<Node> nodes_;
    std::vector<int> scratch_;
    std::vector<int> across_;
    int root_ = -1;
};

} // namespace areko::remesher::curves
