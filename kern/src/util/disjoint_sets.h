#pragma once

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <span>
#include <vector>

namespace areko::remesher::util {

// The numbers 0 … n-1 as sets that can only grow together. Every number
// carries the smallest member of its set, so a set is named by that member.
// Joining rewrites the whole list: meant for few sets.
class DisjointSets {
public:
    explicit DisjointSets(std::size_t count) : smallest_(count) { std::iota(smallest_.begin(), smallest_.end(), 0); }

    // Smallest member of the set of i.
    [[nodiscard]] int smallest(int i) const { return smallest_[static_cast<std::size_t>(i)]; }
    [[nodiscard]] bool joined(int a, int b) const { return smallest(a) == smallest(b); }

    // Makes one set of the sets of a and b; two members of one set change nothing.
    void join(int a, int b)
    {
        const int low = std::min(smallest(a), smallest(b));
        const int high = std::max(smallest(a), smallest(b));
        if (low != high) {
            std::replace(smallest_.begin(), smallest_.end(), high, low);
        }
    }

    // Makes one set of the sets of all members of `group`.
    void joinAll(std::span<const int> group)
    {
        for (std::size_t i = 1; i < group.size(); ++i) {
            join(group[0], group[i]);
        }
    }

private:
    std::vector<int> smallest_;
};

} // namespace areko::remesher::util
