#include "prep_cache.h"

#include <bit>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace areko::remesher::run {

namespace {

// Two independent 64-bit lanes over a stream of words; every word is mixed
// with its position, so order and length count.
class Checksum {
public:
    void add(std::uint64_t word) noexcept
    {
        ++count_;
        lanes_[0] = std::rotl(lanes_[0] ^ mix(word + count_ * 0x9E3779B97F4A7C15ULL), 29) * 0xBF58476D1CE4E5B9ULL;
        lanes_[1] = std::rotl(lanes_[1] ^ mix(word ^ (count_ * 0xD6E8FEB86659FD93ULL)), 31) * 0x94D049BB133111EBULL;
    }

    template <class T>
    void add(std::span<const T> values) noexcept
    {
        add(static_cast<std::uint64_t>(values.size()));
        for (const T value : values) {
            if constexpr (std::is_same_v<T, float>) {
                add(std::bit_cast<std::uint32_t>(value));
            } else {
                add(static_cast<std::uint64_t>(static_cast<std::uint32_t>(value)));
            }
        }
    }

    // An optional part: whether it is there, then its values.
    template <class T>
    void add(const std::optional<std::vector<T>>& values) noexcept
    {
        add(values.has_value() ? 1ULL : 0ULL);
        if (values) {
            add(std::span<const T>(*values));
        }
    }

    void add(const std::string& text) noexcept
    {
        add(static_cast<std::uint64_t>(text.size()));
        for (const char c : text) {
            add(static_cast<std::uint64_t>(static_cast<unsigned char>(c)));
        }
    }

    [[nodiscard]] std::array<std::uint64_t, 2> value() const noexcept { return lanes_; }

private:
    static std::uint64_t mix(std::uint64_t z) noexcept
    {
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

    std::uint64_t count_ = 0;
    std::array<std::uint64_t, 2> lanes_{0x243F6A8885A308D3ULL, 0x13198A2E03707344ULL};
};

std::uint32_t flag(std::uint8_t value) noexcept
{
    return value != 0 ? 1U : 0U;
}

} // namespace

PrepKey PrepKey::of(const input::MeshInput& input, const arm_settings& block, const settings::DerivedSettings& derived)
{
    // InputBuilder::build reads the corner colours, color_density and
    // keep_colors only when the input has corner colours.
    const bool colorsRead = input.cornerColors().has_value();
    Checksum sum;
    sum.add(std::span<const float>(input.points()));
    sum.add(std::span<const std::int32_t>(input.faceSizes()));
    sum.add(std::span<const std::int32_t>(input.cornerPoints()));
    for (const arm_face_values kind : {ARM_FACE_GROUP, ARM_FACE_MATERIAL, ARM_FACE_SMOOTHING_GROUPS}) {
        sum.add(input.faceValues(kind));
    }
    sum.add(input.sharpEdges());
    sum.add(input.cornerNormals());
    sum.add(input.cornerNormalsShared() ? 1ULL : 0ULL);
    sum.add(colorsRead ? 1ULL : 0ULL);
    if (colorsRead) {
        sum.add(input.cornerColors());
    }
    sum.add(input.pointDensity());
    sum.add(static_cast<std::uint64_t>(input.uvSets().size()));
    for (const input::MeshInput::Uv& uv : input.uvSets()) {
        sum.add(uv.name);
        sum.add(std::span<const float>(uv.values));
    }

    PrepKey key;
    key.input = sum.value();
    key.settings = {
        // InputBuilder::build (addColors)
        colorsRead ? 1U | flag(block.color_density) << 1 | flag(block.keep_colors) << 2 : 0U,
        // ImportAttributes::apply
        flag(block.follow_materials),
        flag(block.follow_face_groups),
        flag(block.follow_smoothing_groups),
        flag(block.follow_sharp_edges),
        flag(block.merge_cad_groups),
        flag(block.follow_normal_splits),
        block.weld_seams != 0 ? 1U | flag(block.weld_as_crease) << 1 : 0U,
        // Preparation::run
        flag(block.freeze_borders),
        flag(block.detect_creases),
        flag(block.mirror_at_center),
        // Preparation::run, derived values
        std::bit_cast<std::uint32_t>(derived.quadCount),
        std::bit_cast<std::uint32_t>(derived.edgeLength),
        std::bit_cast<std::uint32_t>(derived.creaseAngleRadians),
        std::bit_cast<std::uint32_t>(derived.denoiseStrength),
    };
    const settings::SymmetrySwitches& s = derived.symmetry;
    key.symmetry = {s.on, s.detect, s.mirrorAtEnd, s.halfMode, s.axis[0], s.axis[1], s.axis[2], s.curvatureSymmetry};
    return key;
}

prep::Prepared PrepCache::copyOf(const prep::Prepared& prepared)
{
    return prep::Prepared{mesh::Mesh(prepared.mesh, true), prepared.groups, prepared.state};
}

std::optional<prep::Prepared> PrepCache::find(const PrepKey& key) const
{
    std::shared_ptr<const Entry> entry;
    {
        const std::lock_guard lock(mutex_);
        entry = entry_;
    }
    if (entry == nullptr || !(entry->key == key)) {
        return std::nullopt;
    }
    return copyOf(entry->prepared);
}

void PrepCache::store(const PrepKey& key, const prep::Prepared& prepared)
{
    auto entry = std::make_shared<const Entry>(Entry{key, copyOf(prepared)});
    const std::lock_guard lock(mutex_);
    entry_ = std::move(entry);
}

} // namespace areko::remesher::run
