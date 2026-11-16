#include "attribute_map.h"

#include <span>

namespace areko::remesher::mesh {

namespace {

// A source counts only when it addresses an existing value; negative or too
// large source numbers are skipped instead of being read.
inline bool usable(int s, std::size_t size) noexcept
{
    return s >= 0 && static_cast<std::size_t>(s) < size;
}

template <class T>
T mixTruncating(const std::vector<T>& values, std::span<const int> sources, const float* weights)
{
    using Traits = ValueTraits<T>;
    T out{};
    const std::size_t n = sources.size();
    if (n == 0) {
        return out;
    }
    std::int64_t sum[Traits::channels] = {};
    if (weights != nullptr) {
        float total = 0.0f;
        for (std::size_t k = 0; k < n; ++k) {
            total += weights[k];
        }
        if (total == 0.0f) {
            return out;
        }
        for (std::size_t k = 0; k < n; ++k) {
            if (!usable(sources[k], values.size())) {
                continue;
            }
            const float share = weights[k] / total;
            for (int ch = 0; ch < Traits::channels; ++ch) {
                sum[ch] += Traits::term(Traits::channel(values[static_cast<std::size_t>(sources[k])], ch) * share);
            }
        }
    } else {
        const float share = 1.0f / static_cast<float>(n);
        for (std::size_t k = 0; k < n; ++k) {
            if (!usable(sources[k], values.size())) {
                continue;
            }
            for (int ch = 0; ch < Traits::channels; ++ch) {
                sum[ch] += Traits::term(Traits::channel(values[static_cast<std::size_t>(sources[k])], ch) * share);
            }
        }
    }
    for (int ch = 0; ch < Traits::channels; ++ch) {
        Traits::setChannel(out, ch, sum[ch]);
    }
    if constexpr (requires { Traits::settle(out, values, sources); }) {
        Traits::settle(out, values, sources);
    }
    return out;
}

// Accumulating rule. An empty row gives the default value of the type; a zero
// weight sum gives ValueTraits<T>::undivided.
template <class T>
T mixAccumulating(const std::vector<T>& values, std::span<const int> sources, const float* weights)
{
    using Traits = ValueTraits<T>;
    const std::size_t n = sources.size();
    if (n == 0) {
        return Traits::resizeDefault();
    }
    typename Traits::Accumulator acc{};
    if (weights != nullptr) {
        float total = 0.0f;
        for (std::size_t k = 0; k < n; ++k) {
            total += weights[k];
            if (usable(sources[k], values.size())) {
                Traits::addWeighted(acc, values[static_cast<std::size_t>(sources[k])], weights[k]);
            }
        }
        return total == 0.0f ? Traits::undivided(acc) : Traits::divide(acc, 1.0f / total);
    }
    for (const int s : sources) {
        if (usable(s, values.size())) {
            Traits::add(acc, values[static_cast<std::size_t>(s)]);
        }
    }
    return Traits::divide(acc, 1.0f / static_cast<float>(n));
}

} // namespace

template <class T>
void TypedMap<T>::interpolateAppend(const Csr& sources, const std::vector<float>* weights)
{
    const int rows = sources.rows();
    std::vector<T> added;
    added.reserve(static_cast<std::size_t>(rows > 0 ? rows : 0));
    for (int r = 0; r < rows; ++r) {
        const auto row = sources.row(r);
        const float* w = weights != nullptr ? weights->data() + sources.starts[static_cast<std::size_t>(r)] : nullptr;
        if constexpr (ValueTraits<T>::truncating) {
            added.push_back(mixTruncating(values_, row, w));
        } else {
            added.push_back(mixAccumulating(values_, row, w));
        }
    }
    values_.insert(values_.end(), added.begin(), added.end());
}

template <class T>
bool TypedMap<T>::reorder(const IndexMap& map)
{
    if (map.newToOld) {
        std::vector<T> out(map.newToOld->size(), ValueTraits<T>::reorderDefault());
        for (std::size_t i = 0; i < out.size(); ++i) {
            const int o = (*map.newToOld)[i];
            if (usable(o, values_.size())) {
                out[i] = values_[static_cast<std::size_t>(o)];
            }
        }
        values_.swap(out);
        return true;
    }
    if (map.oldToNew) {
        std::vector<T> out(static_cast<std::size_t>(map.count), ValueTraits<T>::reorderDefault());
        const std::size_t n = std::min(values_.size(), map.oldToNew->size());
        for (std::size_t i = 0; i < n; ++i) {
            const int j = (*map.oldToNew)[i];
            if (usable(j, out.size())) {
                out[static_cast<std::size_t>(j)] = values_[i];
            }
        }
        values_.swap(out);
        return true;
    }
    return false;
}

template <class T>
bool TypedMap<T>::appendFrom(const AttributeMap& other)
{
    const auto* same = dynamic_cast<const TypedMap<T>*>(&other);
    if (same == nullptr || same->kind() != kind()) {
        return false;
    }
    values_.insert(values_.end(), same->values_.begin(), same->values_.end());
    return true;
}

void BoolMap::interpolateAppend(const Csr& sources, const std::vector<float>* weights)
{
    const std::size_t first = bits_.size();
    const int rows = sources.rows();
    bits_.resize(first + static_cast<std::size_t>(rows > 0 ? rows : 0), false);
    // Sources address the bits that existed before this call.
    const std::size_t size = first;
    for (int r = 0; r < rows; ++r) {
        const auto row = sources.row(r);
        bool on = false;
        if (weights == nullptr) {
            // One or two sources, every usable one set.
            if (!row.empty() && row.size() <= 2) {
                int seen = 0;
                bool all = true;
                for (const int s : row) {
                    if (!usable(s, size)) {
                        continue;
                    }
                    ++seen;
                    all = all && bits_.test(static_cast<std::size_t>(s));
                }
                on = all && seen > 0;
            }
        } else {
            // One or two set sources with a positive weight.
            const float* w = weights->data() + sources.starts[static_cast<std::size_t>(r)];
            int count = 0;
            for (std::size_t k = 0; k < row.size(); ++k) {
                if (usable(row[k], size) && bits_.test(static_cast<std::size_t>(row[k])) && w[k] > 0.0f) {
                    ++count;
                }
            }
            on = count == 1 || count == 2;
        }
        if (on) {
            bits_.set(first + static_cast<std::size_t>(r));
        }
    }
}

bool BoolMap::reorder(const IndexMap& map)
{
    if (map.newToOld) {
        BitArray out(map.newToOld->size());
        for (std::size_t i = 0; i < out.size(); ++i) {
            const int o = (*map.newToOld)[i];
            if (usable(o, bits_.size()) && bits_.test(static_cast<std::size_t>(o))) {
                out.set(i);
            }
        }
        bits_ = std::move(out);
        return true;
    }
    if (map.oldToNew) {
        BitArray out(static_cast<std::size_t>(map.count));
        const std::size_t n = std::min(bits_.size(), map.oldToNew->size());
        for (std::size_t i = 0; i < n; ++i) {
            const int j = (*map.oldToNew)[i];
            if (bits_.test(i) && usable(j, out.size())) {
                out.set(static_cast<std::size_t>(j));
            }
        }
        bits_ = std::move(out);
        return true;
    }
    return false;
}

bool BoolMap::appendFrom(const AttributeMap& other)
{
    const auto* same = dynamic_cast<const BoolMap*>(&other);
    if (same == nullptr || same->kind() != kind()) {
        return false;
    }
    const std::size_t first = bits_.size();
    bits_.resize(first + same->bits_.size());
    for (std::size_t i = 0; i < same->bits_.size(); ++i) {
        if (same->bits_.test(i)) {
            bits_.set(first + i);
        }
    }
    return true;
}

template class TypedMap<int>;
template class TypedMap<std::uint32_t>;
template class TypedMap<std::uint16_t>;
template class TypedMap<std::uint8_t>;
template class TypedMap<float>;
template class TypedMap<double>;
template class TypedMap<Vec3f>;
template class TypedMap<Color3ub>;
template class TypedMap<Color4ub>;
template class TypedMap<FreezeMark>;

} // namespace areko::remesher::mesh
