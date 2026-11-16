#include "ngon_optimizer.h"

#include "edge_dissolve.h"
#include "ngon_finish.h"
#include "ngon_parity_balance.h"
#include "ngon_split_round.h"

#include "mesh/bit_array.h"

#include <bit>
#include <cstdint>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

constexpr int parityRounds = 11;

bool splittable(const Mesh& mesh, int f)
{
    return mesh.topology().mode() == FaceMode::Polygons && mesh.topology().faceSize(f) > 4;
}

} // namespace

NgonSettings NgonSettings::forPart(bool frozenBorderData, bool followBorders)
{
    NgonSettings s;
    if (frozenBorderData) {
        s.balanceParity = true;
        s.smoothNewPoints = true;
        s.smoothingIterations = 4;
        s.smoothingLambda = std::bit_cast<float>(0x3f333333u); // 0.7
    } else if (followBorders) {
        s.smoothNewPoints = true;
    }
    return s;
}

void NgonOptimizer::run(Mesh& mesh, const NgonSettings& settings)
{
    const int startPoints = mesh.pointCount();
    {
        SplitRound round(mesh);
        const int faceCount = mesh.faceCount();
        if (!settings.balanceParity) {
            for (int f = 0; f < faceCount; ++f) {
                if (splittable(mesh, f)) {
                    round.splitFace(f);
                }
            }
            if (round.pendingCount() >= 1) {
                round.commit();
            }
        } else {
            // Even faces first, then rounds of parity pairing and splitting.
            BitArray faces(static_cast<std::size_t>(faceCount), true);
            for (int f = 0; f < faceCount; ++f) {
                if (splittable(mesh, f) && mesh.topology().faceSize(f) % 2 == 0) {
                    round.splitFace(f);
                }
            }
            if (round.pendingCount() > 0) {
                round.commit();
                faces = BitArray(static_cast<std::size_t>(mesh.faceCount()), true);
                round.reset();
            }
            for (int left = parityRounds; left > 0;) {
                ParityBalance::run(round, faces);
                faces.forEachSet([&](std::size_t f) {
                    const int face = static_cast<int>(f);
                    if (!round.done(face) && splittable(mesh, face)) {
                        round.splitFace(face);
                    }
                });
                const int added = round.pendingCount();
                if (added < 1) {
                    break;
                }
                round.commit();
                if (--left == 0) {
                    break;
                }
                const int n = mesh.faceCount();
                faces = BitArray(static_cast<std::size_t>(n));
                for (int f = n - added; f < n; ++f) {
                    faces.set(static_cast<std::size_t>(f));
                }
                round.reset();
            }
        }
    }
    const std::vector<int> doublets = DoubletEdges::find(mesh);
    if (!doublets.empty()) {
        EdgeDissolve::run(mesh, doublets, false);
    }
    if (settings.smoothNewPoints) {
        NewPointSmoother::run(mesh, startPoints, settings.smoothingIterations, settings.smoothingLambda);
    }
}

} // namespace areko::remesher::extract
