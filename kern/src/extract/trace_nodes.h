#pragma once

#include "trace_node.h"

#include <array>
#include <vector>

namespace areko::remesher::extract {

// The nodes of a trace: every node in one list, one list head per triangle
// (newest node first), the counter of output numbers and the stack of nodes
// whose faces are still to be walked. Nodes are named by their index.
class TraceNodes {
public:
    explicit TraceNodes(int triangleCount);

    // A new node in the list of `triangle` (its triangle as well), kind
    // free, without number, flags or links. With `stacked` it is also put on
    // the stack.
    int create(int triangle, mesh::Vec3f position, std::array<float, 2> uv, bool stacked);

    [[nodiscard]] TraceNode& operator[](int node) { return nodes_[static_cast<std::size_t>(node)]; }
    [[nodiscard]] const TraceNode& operator[](int node) const { return nodes_[static_cast<std::size_t>(node)]; }
    [[nodiscard]] int count() const noexcept { return static_cast<int>(nodes_.size()); }
    [[nodiscard]] int triangleCount() const noexcept { return static_cast<int>(heads_.size()); }

    // Calls visit(node) for the nodes in the list of `triangle`, newest first;
    // stops when visit returns true and gives that node, else -1.
    template <class Visit>
    int findIn(int triangle, Visit visit) const
    {
        for (int n = heads_[static_cast<std::size_t>(triangle)]; n >= 0; n = nodes_[static_cast<std::size_t>(n)].older) {
            if (visit(n)) {
                return n;
            }
        }
        return -1;
    }

    // The next output number.
    int takeId() noexcept { return ids_++; }
    [[nodiscard]] int idCount() const noexcept { return ids_; }

    // Stack of nodes (last in, first out).
    [[nodiscard]] bool stackEmpty() const noexcept { return stack_.empty(); }
    int pop();
    // Takes the node put on the stack last off it again.
    void dropTop() { stack_.pop_back(); }

    // State to go back to: the nodes and output numbers so far.
    struct Mark {
        int nodes = 0;
        int ids = 0;
    };
    [[nodiscard]] Mark mark() const noexcept { return {count(), ids_}; }
    // Removes every node made since `since` from the lists of the
    // triangles, from the stack and from the node list; the counter of
    // output numbers goes back to where it stood.
    void discardSince(Mark since);

private:
    std::vector<TraceNode> nodes_;
    std::vector<int> heads_;
    std::vector<int> stack_;
    int ids_ = 0;
};

} // namespace areko::remesher::extract
