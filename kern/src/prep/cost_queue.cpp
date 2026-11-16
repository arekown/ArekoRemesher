#include "cost_queue.h"

#include "mesh/vec3.h"

#include <bit>
#include <cmath>
#include <limits>

namespace areko::remesher::prep {

namespace {

constexpr float ln2 = std::bit_cast<float>(0x3f317218u);
constexpr float log10e = std::bit_cast<float>(0x3ede5bd9u);
constexpr double ln2Divisor = 0.6931471824645996;
constexpr int histogramSize = 256;
constexpr int histogramLimit = 15;

} // namespace

using mesh::truncToInt;

FastLog2::FastLog2()
{
    for (std::uint32_t k = 0; k < table_.size(); ++k) {
        const float f = std::bit_cast<float>(0x3f800000u + (k << 15));
        table_[k] = static_cast<float>(std::log(static_cast<double>(f)) / ln2Divisor);
    }
}

float FastLog2::operator()(float x) const noexcept
{
    const auto bits = std::bit_cast<std::uint32_t>(x);
    return static_cast<float>(static_cast<int>((bits >> 23) & 0xffu) - 127) + mantissaTerm(bits);
}

CostQueue::CostQueue(const std::vector<Collapse>& collapses, float diagonal) : tree_(Order{&collapses})
{
    const int vertices = static_cast<int>(collapses.size());
    binCount_ = vertices / 10;
    heads_.assign(static_cast<std::size_t>(binCount_) + 1, -1);
    next_.assign(static_cast<std::size_t>(vertices), -1);
    prev_ = next_;
    where_.assign(static_cast<std::size_t>(vertices), Where::None);
    costNorm_ = 1e-10;
    logScale_ = 1.0;
    if (vertices == 0) {
        return;
    }
    std::array<int, histogramSize> histogram{};
    for (const Collapse& c : collapses) {
        if (c.cost == std::numeric_limits<float>::max()) {
            continue;
        }
        const int h = truncToInt((log2_(c.cost) * ln2) * log10e) + 128;
        ++histogram[static_cast<std::size_t>(h < 0 ? 0 : (h > 255 ? 255 : h))];
    }
    double low = -32.0;
    int sum = 0;
    for (int k = 0; k < histogramSize; ++k) {
        sum += histogram[static_cast<std::size_t>(k)];
        if (sum > histogramLimit) {
            low = static_cast<double>(k - 128);
            costNorm_ = std::pow(10.0, low);
            break;
        }
    }
    const double reach = static_cast<double>(diagonal) * 0.01;
    const double high = reach * reach;
    logScale_ = (static_cast<double>(binCount_) / (std::log10(high) - low)) * static_cast<double>(log10e);
    for (int v = 0; v < vertices; ++v) {
        const float cost = collapses[static_cast<std::size_t>(v)].cost;
        if (cost != std::numeric_limits<float>::max()) {
            insert(v, bin(cost));
        }
    }
}

int CostQueue::binOfRatio(double ratio) const noexcept
{
    float z = 0.0f;
    if (!(ratio < 1.0)) {
        z = log2_(static_cast<float>(ratio)) * ln2;
    }
    const int f = truncToInt(static_cast<double>(z) * logScale_);
    return f < binCount_ ? f : binCount_;
}

int CostQueue::bin(float cost) const noexcept
{
    return binOfRatio(static_cast<double>(cost) / costNorm_);
}

int CostQueue::binOfDouble(double cost) const noexcept
{
    return binOfRatio(cost / costNorm_);
}

void CostQueue::insert(int v, int bin)
{
    const auto i = static_cast<std::size_t>(v);
    if (bin > current_) {
        const int head = heads_[static_cast<std::size_t>(bin)];
        next_[i] = head;
        if (head >= 0) {
            prev_[static_cast<std::size_t>(head)] = v;
        }
        prev_[i] = -3 - bin;
        heads_[static_cast<std::size_t>(bin)] = v;
        where_[i] = Where::List;
    } else {
        tree_.insert(v);
        where_[i] = Where::Tree;
    }
}

void CostQueue::remove(int v)
{
    const auto i = static_cast<std::size_t>(v);
    if (where_[i] == Where::Tree) {
        tree_.erase(v);
    } else if (where_[i] == Where::List) {
        const int p = prev_[i];
        const int n = next_[i];
        if (p <= -3) {
            heads_[static_cast<std::size_t>(-3 - p)] = n;
        } else {
            next_[static_cast<std::size_t>(p)] = n;
        }
        if (n >= 0) {
            prev_[static_cast<std::size_t>(n)] = p;
        }
    }
    next_[i] = -1;
    prev_[i] = -1;
    where_[i] = Where::None;
}

int CostQueue::front()
{
    while (tree_.empty()) {
        int f = current_ + 1;
        while (f <= binCount_ && heads_[static_cast<std::size_t>(f)] < 0) {
            ++f;
        }
        if (f > binCount_) {
            current_ = binCount_;
            return -1;
        }
        current_ = f;
        for (int v = heads_[static_cast<std::size_t>(f)]; v >= 0;) {
            const int n = next_[static_cast<std::size_t>(v)];
            next_[static_cast<std::size_t>(v)] = -1;
            prev_[static_cast<std::size_t>(v)] = -1;
            tree_.insert(v);
            where_[static_cast<std::size_t>(v)] = Where::Tree;
            v = n;
        }
        heads_[static_cast<std::size_t>(f)] = -1;
    }
    return *tree_.begin();
}

} // namespace areko::remesher::prep
