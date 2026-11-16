#include "closest_face.h"

#include <algorithm>
#include <cfloat>
#include <cmath>

namespace areko::remesher::curves {

using namespace areko::remesher::mesh;

namespace {

// Triangles a node may hold; a split leaving fewer than minSide on one side
// makes a node of at most badSplitLeaf triangles a leaf.
constexpr int leafSize = 5;
constexpr int minSide = 2;
constexpr int badSplitLeaf = 10;

// Nearest point of the segment a-b to p and its distance.
float onSegment(Vec3f p, Vec3f a, Vec3f b, Vec3f& out)
{
    out = nearestOnSegment(p, a, b).point;
    return distance(out, p);
}

// Grows a box (min x, y, z, max x, y, z) by a point, one axis after the other.
void extend(std::array<float, 6>& box, const float* q)
{
    for (int a = 0; a < 3; ++a) {
        const float v = q[a];
        if (box[static_cast<std::size_t>(a)] <= v) {
            if (box[static_cast<std::size_t>(a + 3)] < v) {
                box[static_cast<std::size_t>(a + 3)] = v;
            }
        } else {
            box[static_cast<std::size_t>(a)] = v;
        }
    }
}

// The point lies in the box grown by eps on every side.
bool within(const std::array<float, 6>& box, Vec3f p, float eps)
{
    return box[0] - eps <= p.x && p.x <= box[3] + eps && box[1] - eps <= p.y && p.y <= box[4] + eps &&
           box[2] - eps <= p.z && p.z <= eps + box[5];
}

} // namespace

ClosestFaceSearch::ClosestFaceSearch(const Mesh& mesh) : mesh_(mesh)
{
    const int faces = mesh.faceCount();
    if (faces == 0) {
        return;
    }
    const auto& points = mesh.points();
    const Topology& topology = mesh.topology();
    faceBoxes_.resize(static_cast<std::size_t>(faces));
    order_.resize(static_cast<std::size_t>(faces));
    for (int f = 0; f < faces; ++f) {
        const Vec3f p0 = points[static_cast<std::size_t>(topology.vertex(3 * f))];
        Box box{p0.x, p0.y, p0.z, p0.x, p0.y, p0.z};
        for (int c = 3 * f + 1; c < 3 * f + 3; ++c) {
            const Vec3f q = points[static_cast<std::size_t>(topology.vertex(c))];
            const float xyz[3] = {q.x, q.y, q.z};
            extend(box, xyz);
        }
        faceBoxes_[static_cast<std::size_t>(f)] = box;
        order_[static_cast<std::size_t>(f)] = f;
    }
    const Bounds& b = mesh.bounds();
    nodes_.reserve(static_cast<std::size_t>(faces));
    root_ = build({b.min.x, b.min.y, b.min.z, b.max.x, b.max.y, b.max.z}, 0, faces);
    scratch_ = {};
    across_ = {};
}

ClosestFaceSearch::Box ClosestFaceSearch::unionOf(int begin, int end) const
{
    Box box = faceBoxes_[static_cast<std::size_t>(order_[static_cast<std::size_t>(begin)])];
    for (int i = begin + 1; i < end; ++i) {
        const Box& q = faceBoxes_[static_cast<std::size_t>(order_[static_cast<std::size_t>(i)])];
        extend(box, q.data());
        extend(box, q.data() + 3);
    }
    return box;
}

int ClosestFaceSearch::build(const Box& box, int begin, int end)
{
    const int index = static_cast<int>(nodes_.size());
    nodes_.push_back(Node{box});
    const int count = end - begin;
    if (count <= leafSize) {
        nodes_[static_cast<std::size_t>(index)].leaf = true;
        nodes_[static_cast<std::size_t>(index)].begin = begin;
        nodes_[static_cast<std::size_t>(index)].end = end;
        return index;
    }
    // Longest side: y unless x is longer, z when it is longer than the winner.
    const float ex = box[3] - box[0];
    const float ey = box[4] - box[1];
    const float ez = box[5] - box[2];
    int axis = 1;
    if (ex > ey) {
        axis = ez > ex ? 2 : 0;
    } else if (ez > ey) {
        axis = 2;
    }
    const auto a = static_cast<std::size_t>(axis);
    const float middle = (box[a] + box[a + 3]) * 0.5f;
    scratch_.assign(static_cast<std::size_t>(count), 0);
    across_.clear();
    int left = 0;
    int right = count;
    for (int i = begin; i < end; ++i) {
        const int f = order_[static_cast<std::size_t>(i)];
        const Box& q = faceBoxes_[static_cast<std::size_t>(f)];
        if (!(q[a] < middle)) {
            scratch_[static_cast<std::size_t>(--right)] = f;
        } else if (middle < q[a + 3]) {
            across_.push_back(f);
        } else {
            scratch_[static_cast<std::size_t>(left++)] = f;
        }
    }
    for (const int f : across_) {
        const Box& q = faceBoxes_[static_cast<std::size_t>(f)];
        if (middle - q[a] >= q[a + 3] - middle) {
            scratch_[static_cast<std::size_t>(left++)] = f;
        } else {
            scratch_[static_cast<std::size_t>(--right)] = f;
        }
    }
    std::copy(scratch_.begin(), scratch_.end(), order_.begin() + begin);
    int split = begin + left;
    const int middleIndex = static_cast<int>((static_cast<unsigned>(end) + static_cast<unsigned>(begin)) >> 1);
    if (split < begin + minSide) {
        split = middleIndex;
        if (count <= badSplitLeaf) {
            nodes_[static_cast<std::size_t>(index)].leaf = true;
            nodes_[static_cast<std::size_t>(index)].begin = begin;
            nodes_[static_cast<std::size_t>(index)].end = end;
            return index;
        }
    }
    if (end - minSide < split) {
        split = middleIndex;
        if (count <= badSplitLeaf) {
            nodes_[static_cast<std::size_t>(index)].leaf = true;
            nodes_[static_cast<std::size_t>(index)].begin = begin;
            nodes_[static_cast<std::size_t>(index)].end = end;
            return index;
        }
    }
    const Box leftBox = unionOf(begin, split);
    const Box rightBox = unionOf(split, end);
    const int l = build(leftBox, begin, split);
    const int r = build(rightBox, split, end);
    nodes_[static_cast<std::size_t>(index)].left = l;
    nodes_[static_cast<std::size_t>(index)].right = r;
    return index;
}

ClosestFaceSearch::Hit ClosestFaceSearch::onTriangle(int f, Vec3f p) const
{
    const auto& points = mesh_.points();
    const Topology& topology = mesh_.topology();
    const Vec3f p0 = points[static_cast<std::size_t>(topology.vertex(3 * f))];
    const Vec3f p1 = points[static_cast<std::size_t>(topology.vertex(3 * f + 1))];
    const Vec3f p2 = points[static_cast<std::size_t>(topology.vertex(3 * f + 2))];
    Hit hit;
    hit.face = f;
    Vec3f n = mesh_.triangleGeometry().normals()[static_cast<std::size_t>(f)];
    const float nl = length(n);
    if (!(nl != 0.0f)) {
        hit.distance = FLT_MAX;
        return hit;
    }
    n = {n.x / nl, n.y / nl, n.z / nl};
    // Drop p onto the plane of the triangle.
    const float h = dot(n, p) - dot(n, p0);
    const Vec3f q = p - n * h;
    const Vec3f e1 = p1 - p0;
    const Vec3f e2 = p2 - p0;
    const Vec3f w = q - p0;
    const float d11 = dot(e1, e1);
    const float d12 = dot(e1, e2);
    const float d22 = dot(e2, e2);
    const float b1 = dot(e1, w);
    const float b2 = dot(e2, w);
    const float inv = 1.0f / (d11 * d22 - d12 * d12);
    const float u = (d11 * b2 - b1 * d12) * inv; // along e2
    const float v = (d22 * b1 - d12 * b2) * inv; // along e1
    if (v < 0.0f || u < 0.0f || 1.0f < v + u) {
        Vec3f s1;
        Vec3f s2;
        Vec3f s3;
        const float d1 = onSegment(p, p0, p1, s1);
        const float d2 = onSegment(p, p1, p2, s2);
        const float d3 = onSegment(p, p2, p0, s3);
        if (d2 <= d1) {
            hit.point = d3 <= d2 ? s3 : s2;
        } else {
            hit.point = d3 <= d1 ? s3 : s1;
        }
    } else {
        hit.point = q;
    }
    hit.distance = distance(hit.point, p);
    return hit;
}

ClosestFaceSearch::Hit ClosestFaceSearch::nearest(Vec3f p, float radius) const
{
    Hit best;
    best.distance = radius;
    if (root_ == -1) {
        return best;
    }
    std::vector<int> stack{root_};
    while (!stack.empty()) {
        const Node& node = nodes_[static_cast<std::size_t>(stack.back())];
        stack.pop_back();
        if (!within(node.box, p, best.distance)) {
            continue;
        }
        if (!node.leaf) {
            // The right child comes off the stack first.
            stack.push_back(node.left);
            stack.push_back(node.right);
            continue;
        }
        for (int i = node.begin; i < node.end; ++i) {
            const Hit hit = onTriangle(order_[static_cast<std::size_t>(i)], p);
            if (hit.distance < best.distance) {
                best = hit;
            }
        }
    }
    return best;
}

} // namespace areko::remesher::curves
