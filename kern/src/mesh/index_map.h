#pragma once

#include <cstdint>
#include <optional>
#include <vector>

namespace areko::remesher::mesh {

// Renumbering of one element kind after a mesh change. -1 means "no element".
// When both directions are given, new -> old is preferred.
struct IndexMap {
    std::optional<std::vector<int>> oldToNew;
    std::optional<std::vector<int>> newToOld;
    int count = 0;           // number of new elements (used with oldToNew)
    std::uint8_t state = 0;  // 1 = element kind unchanged; other values are not told apart

    [[nodiscard]] bool unchanged() const noexcept { return (state & 3) == 1; }
    [[nodiscard]] bool present() const noexcept { return oldToNew.has_value() || newToOld.has_value(); }
};

} // namespace areko::remesher::mesh
