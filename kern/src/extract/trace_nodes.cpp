#include "trace_nodes.h"

#include <algorithm>

namespace areko::remesher::extract {

TraceNodes::TraceNodes(int triangleCount) : heads_(static_cast<std::size_t>(triangleCount), -1) {}

int TraceNodes::create(int triangle, mesh::Vec3f position, std::array<float, 2> uv, bool stacked)
{
    const int node = count();
    TraceNode& n = nodes_.emplace_back();
    n.position = position;
    n.triangle = triangle;
    n.uv = uv;
    int& head = heads_[static_cast<std::size_t>(triangle)];
    n.older = head;
    head = node;
    if (stacked) {
        stack_.push_back(node);
    }
    return node;
}

int TraceNodes::pop()
{
    const int node = stack_.back();
    stack_.pop_back();
    return node;
}

void TraceNodes::discardSince(Mark since)
{
    // Newer nodes stand in front of older ones in every list.
    for (int& head : heads_) {
        while (head >= since.nodes) {
            head = nodes_[static_cast<std::size_t>(head)].older;
        }
    }
    std::erase_if(stack_, [&](int node) { return node >= since.nodes; });
    nodes_.resize(static_cast<std::size_t>(since.nodes));
    ids_ = since.ids;
}

} // namespace areko::remesher::extract
