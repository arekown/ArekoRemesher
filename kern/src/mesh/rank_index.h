#pragma once

#include <algorithm>
#include <numeric>
#include <vector>

namespace areko::remesher::mesh {

// Numbers equal keys: element i gets the rank of its key among the distinct
// keys in ascending order. `representative[r]` is the first element (in sorted
// order, ties by element index) carrying rank r.
template <class Key>
std::vector<int> rankByKey(const std::vector<Key>& keys, std::vector<int>& representative)
{
    std::vector<int> order(keys.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&keys](int a, int b) {
        return keys[a] < keys[b] || (!(keys[b] < keys[a]) && a < b);
    });
    std::vector<int> rank(keys.size(), 0);
    representative.clear();
    for (std::size_t i = 0; i < order.size(); ++i) {
        if (i == 0 || keys[order[i - 1]] < keys[order[i]]) {
            representative.push_back(order[i]);
        }
        rank[order[i]] = static_cast<int>(representative.size()) - 1;
    }
    return rank;
}

} // namespace areko::remesher::mesh
