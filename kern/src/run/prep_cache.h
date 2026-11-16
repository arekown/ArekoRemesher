#pragma once

#include "areko_remesher.h"
#include "input/mesh_input.h"
#include "prep/preparation.h"
#include "settings/settings.h"

#include <array>
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>

namespace areko::remesher::run {

// What the stages "input" and "prepared" read: a checksum of the input parts
// they read and every setting they read, as it reaches them.
//  - input: points, faces, UV sets, normals, face values, sharp edges, density
//    and, when colours are built, the corner colours (not the guide curves);
//  - settings block: color_density and keep_colors (only with corner colours),
//    the follow_* switches of the import, merge_cad_groups, weld_seams,
//    weld_as_crease (only with weld_seams), freeze_borders, detect_creases,
//    mirror_at_center;
//  - derived values: quad count, edge length, crease angle, denoise strength
//    and the symmetry switches.
struct PrepKey {
    std::array<std::uint64_t, 2> input{};
    std::array<std::uint32_t, 15> settings{};
    std::array<std::uint8_t, 8> symmetry{};

    friend bool operator==(const PrepKey&, const PrepKey&) = default;

    [[nodiscard]] static PrepKey of(const input::MeshInput& input, const arm_settings& block,
                                    const settings::DerivedSettings& derived);
};

// The remembered preparation of one key: the state after the stage
// "prepared", before the guide curves are looked at. Shared by the contexts of
// one host session (any thread); every call copies, so a run never changes it.
class PrepCache {
public:
    // A copy of the remembered preparation when the key matches.
    [[nodiscard]] std::optional<prep::Prepared> find(const PrepKey& key) const;
    // Remembers a copy of `prepared` for `key`, in place of the earlier one.
    void store(const PrepKey& key, const prep::Prepared& prepared);

    // A copy that keeps the derived caches of the mesh the copy would rebuild
    // the same way.
    [[nodiscard]] static prep::Prepared copyOf(const prep::Prepared& prepared);

private:
    struct Entry {
        PrepKey key;
        prep::Prepared prepared;
    };

    mutable std::mutex mutex_;
    std::shared_ptr<const Entry> entry_;
};

} // namespace areko::remesher::run
