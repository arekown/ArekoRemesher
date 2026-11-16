#include "seam_rows.h"

#include "mesh/vec3.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace areko::remesher::param {

using namespace areko::remesher::mesh;

namespace {

// Whether the two sides count as lying on one grid line: `minor` and `major`
// are the smaller and the larger part of their offset, `length` is how far
// the chain runs in the coordinate of the smaller part.
[[nodiscard]] bool onOneLine(float minor, float major, float length)
{
    const auto skew = [&] { return wide(minor / (major * 0.75f + length * 0.25f)); };
    if (minor < 1.0f) {
        return (major > 4.0f && length > 6.0f) || (major > 6.0f && length > 4.0f)
               || (major > 2.0f && length > 2.0f && skew() < 0.02);
    }
    return minor < 2.0f && major > 14.0f && length > 10.0f && skew() < 0.08;
}

} // namespace

void SeamRows::add(const SystemContext& context, const UvOutput::Coordinates& uv, SoftSystem& system)
{
    const Chains& chains = context.chains();
    const BitArray& diskCuts = context.graph().diskCutBits;
    for (int chain = 0; chain < chains.count(); ++chain) {
        const auto i = static_cast<std::size_t>(chain);
        if ((truncToInt(chains.rotation[i]) & 3) != 0) {
            continue;
        }
        const std::vector<int>& sideA = chains.sideA[i];
        if (sideA.empty()) {
            continue;
        }
        if (diskCuts.testBounded(static_cast<std::size_t>(chains.edges.row(chain)[0]))) {
            continue;
        }
        const int firstA = sideA.front();
        const int firstB = chains.sideB[i].front();
        const std::array<float, 2>& a = uv[static_cast<std::size_t>(firstA)];
        const std::array<float, 2>& b = uv[static_cast<std::size_t>(firstB)];
        const std::array<float, 2>& end = uv[static_cast<std::size_t>(sideA.back())];
        const float offsetU = std::fabs(b[0] - a[0]);
        const float offsetV = std::fabs(b[1] - a[1]);
        // The coordinate of the smaller part; the second one unless the
        // offset in it is the larger one.
        const std::size_t m = offsetV > offsetU ? 0 : 1;
        const float minor = m == 1 ? offsetV : offsetU;
        const float major = m == 1 ? offsetU : offsetV;
        const float length = std::fabs(end[m] - a[m]);
        if (!onOneLine(minor, major, length)) {
            continue;
        }
        system.addHard(TermList{{2 * firstA + static_cast<int>(m), 1.0}, {2 * firstB + static_cast<int>(m), -1.0}},
                       0.0);
    }
}

} // namespace areko::remesher::param
