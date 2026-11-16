#pragma once

#include <array>
#include <cstdint>
#include <set>
#include <vector>

namespace areko::remesher::prep {

// Base-2 logarithm from a table of 256 mantissa steps plus the exponent.
class FastLog2 {
public:
    FastLog2();
    [[nodiscard]] float operator()(float x) const noexcept;

private:
    [[nodiscard]] float mantissaTerm(std::uint32_t bits) const noexcept { return table_[(bits & 0x7fffffu) >> 15]; }

    std::array<float, 256> table_{};
};

// Best collapse of a vertex: target vertex and cost.
struct Collapse {
    int target = -1;
    float cost = 0.0f;
};

// Queue of vertices ordered by (cost, vertex number). Costs are spread over
// logarithmic bins; bins up to the current one live in an ordered tree, the
// others in lists. The front of the tree is always the cheapest vertex. The
// costs are read live from the collapse table, so a vertex must leave the
// queue before its cost changes.
class CostQueue {
public:
    // Sets up the bins from the start costs (FLT_MAX = not queued) and the
    // bounding box diagonal.
    CostQueue(const std::vector<Collapse>& collapses, float diagonal);
    CostQueue(const CostQueue&) = delete;
    CostQueue& operator=(const CostQueue&) = delete;

    [[nodiscard]] int bin(float cost) const noexcept;
    [[nodiscard]] int binOfDouble(double cost) const noexcept;
    [[nodiscard]] int binCount() const noexcept { return binCount_; }
    [[nodiscard]] bool queued(int v) const noexcept { return where_[static_cast<std::size_t>(v)] != Where::None; }

    void insert(int v, int bin);
    void remove(int v);
    // Cheapest vertex (it stays queued), or -1 when the queue is empty.
    [[nodiscard]] int front();

private:
    enum class Where : std::uint8_t { None, List, Tree };
    struct Order {
        const std::vector<Collapse>* collapses;
        bool operator()(int a, int b) const
        {
            const float ca = (*collapses)[static_cast<std::size_t>(a)].cost;
            const float cb = (*collapses)[static_cast<std::size_t>(b)].cost;
            return ca < cb || (ca == cb && a < b);
        }
    };
    [[nodiscard]] int binOfRatio(double ratio) const noexcept;

    FastLog2 log2_;
    int binCount_ = 0;
    double costNorm_ = 1e-10;
    double logScale_ = 1.0;
    int current_ = 0;
    std::vector<int> heads_;
    std::vector<int> next_;
    std::vector<int> prev_;
    std::vector<Where> where_;
    std::set<int, Order> tree_;
};

} // namespace areko::remesher::prep
