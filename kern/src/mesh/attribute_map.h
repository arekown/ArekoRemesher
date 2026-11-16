#pragma once

#include "bit_array.h"
#include "csr.h"
#include "index_map.h"
#include "vec3.h"

#include <memory>
#include <span>
#include <string>
#include <vector>

namespace areko::remesher::mesh {

// Element a map belongs to; the numbers are part of the data model.
enum class ElementKind : int {
    Face = 1,
    Point = 2,
    Corner = 3,
    UvPoint = 4,
    Edge = 5,
    UvEdge = 6,
    HalfEdge = 7
};

// One named attribute with one value per element.
class AttributeMap {
public:
    AttributeMap(std::string name, ElementKind kind) : name_(std::move(name)), kind_(kind) {}
    virtual ~AttributeMap() = default;
    AttributeMap(const AttributeMap&) = default;
    AttributeMap& operator=(const AttributeMap&) = delete;

    [[nodiscard]] const std::string& name() const noexcept { return name_; }
    [[nodiscard]] ElementKind kind() const noexcept { return kind_; }

    [[nodiscard]] virtual std::size_t size() const noexcept = 0;
    // Appends one value per row of `sources`, mixed from the listed values.
    // `weights` (may be null) runs parallel to sources.entries.
    virtual void interpolateAppend(const Csr& sources, const std::vector<float>* weights) = 0;
    // Renumbers the values; false (map unchanged) when `map` has no direction.
    virtual bool reorder(const IndexMap& map) = 0;
    // Appends the values of a map of the same type and kind; false otherwise.
    virtual bool appendFrom(const AttributeMap& other) = 0;
    // Appends default values or cuts to `size`.
    virtual void resize(std::size_t size) = 0;
    [[nodiscard]] virtual std::unique_ptr<AttributeMap> clone() const = 0;

private:
    std::string name_;
    ElementKind kind_;
};

// Rules per value type: how rows are mixed and which value fills gaps.
//  - "truncating" types sum the truncated, weighted terms per channel,
//  - the others accumulate in `Accumulator` (addWeighted / add per usable
//    source) and divide once; `undivided` is the result for a zero weight sum.
template <class T>
struct ValueTraits;

template <>
struct ValueTraits<int> {
    static constexpr bool truncating = false;
    static int reorderDefault() { return 0; }
    static int resizeDefault() { return 0; }
    // Unweighted sums wrap around in 32 bits; weighted sums are truncated after every step.
    using Accumulator = std::uint32_t;
    static void addWeighted(Accumulator& acc, int v, float w)
    {
        acc = static_cast<Accumulator>(
            truncToInt(static_cast<float>(static_cast<std::int32_t>(acc)) + static_cast<float>(v) * w));
    }
    static void add(Accumulator& acc, int v) { acc += static_cast<Accumulator>(v); }
    static int divide(Accumulator acc, float inverse)
    {
        return truncToInt(inverse * static_cast<float>(static_cast<std::int32_t>(acc)));
    }
    static int undivided(Accumulator acc) { return static_cast<std::int32_t>(acc); }
};
template <>
struct ValueTraits<float> {
    static constexpr bool truncating = false;
    static float reorderDefault() { return 0.0f; }
    static float resizeDefault() { return 0.0f; }
    using Accumulator = float;
    static void addWeighted(Accumulator& acc, float v, float w) { acc += v * w; }
    static void add(Accumulator& acc, float v) { acc += v; }
    static float divide(Accumulator acc, float inverse) { return inverse * acc; }
    static float undivided(Accumulator acc) { return acc; }
};
template <>
struct ValueTraits<double> {
    static constexpr bool truncating = false;
    static double reorderDefault() { return 0.0; }
    static double resizeDefault() { return 0.0; }
    using Accumulator = double;
    static void addWeighted(Accumulator& acc, double v, float w) { acc += static_cast<double>(w) * v; }
    static void add(Accumulator& acc, double v) { acc += v; }
    static double divide(Accumulator acc, float inverse) { return static_cast<double>(inverse) * acc; }
    static double undivided(Accumulator acc) { return acc; }
};
template <>
struct ValueTraits<Vec3f> {
    static constexpr bool truncating = false;
    static Vec3f reorderDefault() { return {}; }
    static Vec3f resizeDefault() { return {}; }
    using Accumulator = Vec3f;
    static void addWeighted(Accumulator& acc, const Vec3f& v, float w)
    {
        acc.x += v.x * w;
        acc.y += v.y * w;
        acc.z += v.z * w;
    }
    static void add(Accumulator& acc, const Vec3f& v)
    {
        acc.x += v.x;
        acc.y += v.y;
        acc.z += v.z;
    }
    static Vec3f divide(const Accumulator& acc, float inverse) { return {inverse * acc.x, inverse * acc.y, inverse * acc.z}; }
    // A zero weight sum gives the zero vector.
    static Vec3f undivided(const Accumulator&) { return {}; }
};
template <>
struct ValueTraits<std::uint32_t> {
    static constexpr bool truncating = true;
    static constexpr int channels = 1;
    static std::uint32_t reorderDefault() { return 0; }
    static std::uint32_t resizeDefault() { return 0; }
    static float channel(std::uint32_t v, int) { return static_cast<float>(static_cast<std::int64_t>(v)); }
    static std::int64_t term(float v) { return truncToInt64(v); }
    static void setChannel(std::uint32_t& v, int, std::int64_t s) { v = static_cast<std::uint32_t>(s); }
};
template <>
struct ValueTraits<std::uint16_t> {
    static constexpr bool truncating = true;
    static constexpr int channels = 1;
    static std::uint16_t reorderDefault() { return 0; }
    static std::uint16_t resizeDefault() { return 0; }
    static float channel(std::uint16_t v, int) { return static_cast<float>(static_cast<int>(v)); }
    static std::int64_t term(float v) { return truncToInt(v); }
    static void setChannel(std::uint16_t& v, int, std::int64_t s) { v = static_cast<std::uint16_t>(s); }
};
template <>
struct ValueTraits<std::uint8_t> {
    static constexpr bool truncating = true;
    static constexpr int channels = 1;
    static std::uint8_t reorderDefault() { return 0; }
    static std::uint8_t resizeDefault() { return 0; }
    static float channel(std::uint8_t v, int) { return static_cast<float>(static_cast<int>(v)); }
    static std::int64_t term(float v) { return truncToInt(v); }
    static void setChannel(std::uint8_t& v, int, std::int64_t s) { v = static_cast<std::uint8_t>(s); }
};
template <>
struct ValueTraits<Color3ub> {
    static constexpr bool truncating = true;
    static constexpr int channels = 3;
    static Color3ub reorderDefault() { return {}; }
    static Color3ub resizeDefault() { return {}; }
    static float channel(const Color3ub& v, int ch) { return static_cast<float>((&v.r)[ch]); }
    static std::int64_t term(float v) { return truncToInt(v); }
    static void setChannel(Color3ub& v, int ch, std::int64_t s) { (&v.r)[ch] = static_cast<std::uint8_t>(s); }
};
template <>
struct ValueTraits<Color4ub> {
    static constexpr bool truncating = true;
    static constexpr int channels = 4;
    // Renumbering fills with opaque black, growing fills with all four channels zero.
    static Color4ub reorderDefault() { return {0, 0, 0, 255}; }
    static Color4ub resizeDefault() { return {0, 0, 0, 0}; }
    static float channel(const Color4ub& v, int ch) { return static_cast<float>((&v.r)[ch]); }
    static std::int64_t term(float v) { return truncToInt(v); }
    static void setChannel(Color4ub& v, int ch, std::int64_t s) { (&v.r)[ch] = static_cast<std::uint8_t>(s); }
};

// Freeze mark of a point: 0 = none, 1 = kept border point, 2 = between kept points.
struct FreezeMark {
    std::uint8_t value = 0;
    friend bool operator==(const FreezeMark&, const FreezeMark&) = default;
};
template <>
struct ValueTraits<FreezeMark> {
    static constexpr bool truncating = true;
    static constexpr int channels = 1;
    static FreezeMark reorderDefault() { return {}; }
    static FreezeMark resizeDefault() { return {}; }
    static float channel(const FreezeMark& v, int) { return static_cast<float>(static_cast<int>(v.value)); }
    static std::int64_t term(float v) { return truncToInt(v); }
    static void setChannel(FreezeMark& v, int, std::int64_t s) { v.value = static_cast<std::uint8_t>(s); }
    // A point mixed only from kept border points is kept as well.
    static void settle(FreezeMark& out, const std::vector<FreezeMark>& values, std::span<const int> sources)
    {
        bool kept = !sources.empty();
        for (const int s : sources) {
            kept = kept && s >= 0 && static_cast<std::size_t>(s) < values.size() &&
                   values[static_cast<std::size_t>(s)].value == 1;
        }
        if (kept) {
            out.value = 1;
        }
    }
};

template <class T>
class TypedMap final : public AttributeMap {
public:
    TypedMap(std::string name, ElementKind kind, std::vector<T> values = {})
        : AttributeMap(std::move(name), kind), values_(std::move(values))
    {
    }

    [[nodiscard]] std::vector<T>& values() noexcept { return values_; }
    [[nodiscard]] const std::vector<T>& values() const noexcept { return values_; }

    [[nodiscard]] std::size_t size() const noexcept override { return values_.size(); }
    void interpolateAppend(const Csr& sources, const std::vector<float>* weights) override;
    bool reorder(const IndexMap& map) override;
    bool appendFrom(const AttributeMap& other) override;
    void resize(std::size_t size) override { values_.resize(size, ValueTraits<T>::resizeDefault()); }
    // Resize with an explicit fill value.
    void resizeFilled(std::size_t size, const T& fill) { values_.resize(size, fill); }
    [[nodiscard]] std::unique_ptr<AttributeMap> clone() const override { return std::make_unique<TypedMap<T>>(*this); }

private:
    std::vector<T> values_;
};

// One bit per element.
class BoolMap final : public AttributeMap {
public:
    BoolMap(std::string name, ElementKind kind, BitArray bits = {})
        : AttributeMap(std::move(name), kind), bits_(std::move(bits))
    {
    }

    [[nodiscard]] BitArray& bits() noexcept { return bits_; }
    [[nodiscard]] const BitArray& bits() const noexcept { return bits_; }

    [[nodiscard]] std::size_t size() const noexcept override { return bits_.size(); }
    void interpolateAppend(const Csr& sources, const std::vector<float>* weights) override;
    bool reorder(const IndexMap& map) override;
    bool appendFrom(const AttributeMap& other) override;
    void resize(std::size_t size) override { bits_.resize(size, false); }
    [[nodiscard]] std::unique_ptr<AttributeMap> clone() const override { return std::make_unique<BoolMap>(*this); }

private:
    BitArray bits_;
};

extern template class TypedMap<int>;
extern template class TypedMap<std::uint32_t>;
extern template class TypedMap<std::uint16_t>;
extern template class TypedMap<std::uint8_t>;
extern template class TypedMap<float>;
extern template class TypedMap<double>;
extern template class TypedMap<Vec3f>;
extern template class TypedMap<Color3ub>;
extern template class TypedMap<Color4ub>;
extern template class TypedMap<FreezeMark>;

} // namespace areko::remesher::mesh
