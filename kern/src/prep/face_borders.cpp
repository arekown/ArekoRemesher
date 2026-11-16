#include "face_borders.h"

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

template <class T, class Differs>
BitArray borders(const Adjacency& adjacency, std::span<const T> values, bool withBoundary, Differs differs)
{
    BitArray out(static_cast<std::size_t>(adjacency.edgeCount()));
    auto valueAt = [&](int corner, T& value) {
        const int f = adjacency.cornerFace(corner);
        if (f < 0 || static_cast<std::size_t>(f) >= values.size()) {
            return false;
        }
        value = values[static_cast<std::size_t>(f)];
        return true;
    };
    for (int e = 0; e < adjacency.edgeCount(); ++e) {
        const auto corners = adjacency.edgeCorners(e);
        if (corners.size() == 1) {
            if (withBoundary) {
                out.set(static_cast<std::size_t>(e));
            }
            continue;
        }
        T first{};
        if (corners.empty() || !valueAt(corners[0], first)) {
            continue;
        }
        for (std::size_t k = 1; k < corners.size(); ++k) {
            T other{};
            if (valueAt(corners[k], other) && differs(first, other)) {
                out.set(static_cast<std::size_t>(e));
                break;
            }
        }
    }
    return out;
}

} // namespace

BitArray FaceBorders::valueBorders(const Adjacency& adjacency, std::span<const int> faceValues, bool withBoundary)
{
    return borders(adjacency, faceValues, withBoundary, [](int a, int b) { return a != b; });
}

BitArray FaceBorders::smoothingBorders(const Adjacency& adjacency, std::span<const std::uint32_t> faceMasks)
{
    return borders(adjacency, faceMasks, false, [](std::uint32_t a, std::uint32_t b) { return (a & b) == 0; });
}

} // namespace areko::remesher::prep
