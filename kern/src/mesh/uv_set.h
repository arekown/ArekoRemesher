#pragma once

#include <array>
#include <string>
#include <vector>

namespace areko::remesher::mesh {

// One UV set: every corner names a UV point; UV points carry the coordinates.
struct UvSet {
    std::string name;
    std::vector<int> cornerUv;
    std::vector<std::array<float, 2>> coords;

    void clear()
    {
        cornerUv.clear();
        coords.clear();
    }
};

} // namespace areko::remesher::mesh
