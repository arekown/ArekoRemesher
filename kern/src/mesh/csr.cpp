#include "csr.h"

namespace areko::remesher::mesh {

Csr Csr::fromValues(std::span<const int> values, int rowCount)
{
    Csr csr;
    csr.starts.assign(static_cast<std::size_t>(rowCount) + 1, 0);
    for (const int v : values) {
        if (v >= 0) {
            ++csr.starts[static_cast<std::size_t>(v) + 1];
        }
    }
    for (int r = 0; r < rowCount; ++r) {
        csr.starts[r + 1] += csr.starts[r];
    }
    csr.entries.resize(static_cast<std::size_t>(csr.starts[rowCount]));
    std::vector<int> fill(csr.starts.begin(), csr.starts.end() - 1);
    for (std::size_t i = 0; i < values.size(); ++i) {
        const int v = values[i];
        if (v >= 0) {
            csr.entries[static_cast<std::size_t>(fill[v]++)] = static_cast<int>(i);
        }
    }
    return csr;
}

} // namespace areko::remesher::mesh
