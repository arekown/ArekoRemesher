#include "duplicate_faces.h"

#include "mesh/modification.h"

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

int DuplicateFaces::remove(Mesh& mesh)
{
    const Topology& topology = mesh.topology();
    const Adjacency& adjacency = mesh.adjacency();
    BitArray removed(static_cast<std::size_t>(topology.faceCount()));
    int count = 0;
    auto at = [&](int f, int position) { return topology.vertex(topology.faceBegin(f) + position); };
    for (int v = 0; v < adjacency.pointCount(); ++v) {
        const auto corners = adjacency.pointCorners(v);
        const std::size_t k = corners.size();
        for (std::size_t i = 0; i + 1 < k; ++i) {
            const int fi = adjacency.cornerFace(corners[i]);
            if (removed.test(static_cast<std::size_t>(fi))) {
                continue;
            }
            const int pi = corners[i] - topology.faceBegin(fi);
            const int n = topology.faceSize(fi);
            for (std::size_t j = i + 1; j < k; ++j) {
                const int fj = adjacency.cornerFace(corners[j]);
                if (removed.test(static_cast<std::size_t>(fj)) || topology.faceSize(fj) != n) {
                    continue;
                }
                const int pj = corners[j] - topology.faceBegin(fj);
                bool same = true;
                for (int s = 1; s <= n && same; ++s) {
                    same = at(fi, (pi + s) % n) == at(fj, (pj + s) % n);
                }
                if (!same) {
                    same = true;
                    for (int s = 1; s < n && same; ++s) {
                        same = at(fi, (pi + s) % n) == at(fj, ((pj - s) % n + n) % n);
                    }
                }
                if (same) {
                    removed.set(static_cast<std::size_t>(fj));
                    ++count;
                }
            }
        }
    }
    if (count > 0) {
        Modification change;
        change.deleteFaces = std::move(removed);
        MeshModifier::apply(mesh, change);
    }
    return count;
}

} // namespace areko::remesher::prep
