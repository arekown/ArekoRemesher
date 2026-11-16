#include "border_lock.h"

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

BorderLock::Options BorderLock::optionsFor(const settings::SymmetrySwitches& symmetry)
{
    Options options;
    options.halfMode = symmetry.halfMode;
    const bool anyAxis = symmetry.axis[0] || symmetry.axis[1] || symmetry.axis[2];
    const bool givenHalf = symmetry.halfMode && !symmetry.on && !anyAxis;
    for (int i = 0; i < 3; ++i) {
        options.axis[i] = givenHalf || symmetry.axis[i];
    }
    return options;
}

void BorderLock::assign(parts::Part& part, const Options& options)
{
    part.lockedBorderEdges.clear();
    if (!options.halfMode || !(options.axis[0] || options.axis[1] || options.axis[2])) {
        return;
    }
    const Adjacency& adjacency = part.mesh.adjacency();
    const std::vector<Vec3f>& points = part.mesh.points();
    BitArray locked(static_cast<std::size_t>(adjacency.edgeCount()));
    bool any = false;
    for (const int e : adjacency.boundaryEdges()) {
        const auto ends = adjacency.edgeVertices(e);
        const Vec3f a = points[static_cast<std::size_t>(ends[0])];
        const Vec3f b = points[static_cast<std::size_t>(ends[1])];
        for (int i = 0; i < 3; ++i) {
            if (options.axis[i] && axisValue(a, i) == axisValue(b, i) &&
                axisValue(a, i) == axisValue(options.planeMax, i)) {
                locked.set(static_cast<std::size_t>(e));
                any = true;
            }
        }
    }
    if (any) {
        part.lockedBorderEdges = std::move(locked);
    }
}

} // namespace areko::remesher::field
