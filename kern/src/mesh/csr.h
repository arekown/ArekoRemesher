#pragma once

#include <span>
#include <vector>

namespace areko::remesher::mesh {

// Compressed rows: row r holds entries[starts[r] .. starts[r+1]-1].
// A table counts as "not built" while `starts` has exactly one element.
struct Csr {
    std::vector<int> starts{0};
    std::vector<int> entries;

    [[nodiscard]] bool built() const noexcept { return starts.size() > 1; }
    [[nodiscard]] int rows() const noexcept { return static_cast<int>(starts.size()) - 1; }
    [[nodiscard]] int rowSize(int r) const { return starts[r + 1] - starts[r]; }
    [[nodiscard]] std::span<const int> row(int r) const
    {
        return {entries.data() + starts[r], static_cast<std::size_t>(starts[r + 1] - starts[r])};
    }
    void reset()
    {
        starts.assign(1, 0);
        entries.clear();
    }

    // Row v lists every index i with values[i] == v, ascending (stable counting
    // sort). Negative values are skipped.
    [[nodiscard]] static Csr fromValues(std::span<const int> values, int rowCount);
};

} // namespace areko::remesher::mesh
