#pragma once

#include "mesh/boundary_loops.h"
#include "mesh/mesh.h"

#include <cstddef>
#include <span>
#include <vector>

namespace areko::remesher::prep {

// Marks the border points that stay in place when borders are frozen (point
// map of FreezeMark values, 1 = kept border point). With `onlyIfMissing` a
// complete map that is already there stays as it is.
class FreezeMarks {
public:
    static void mark(mesh::Mesh& mesh, bool onlyIfMissing);
    [[nodiscard]] static const mesh::TypedMap<mesh::FreezeMark>* find(const mesh::Mesh& mesh);
    // Positions of a border loop (its points in walking order) whose points
    // are kept, in walking order starting behind position 0, so position 0
    // comes last. Empty for a loop with fewer than two kept points: such a
    // loop has no chord between kept points.
    [[nodiscard]] static std::vector<int> keptPositions(std::span<const int> loopPoints,
                                                        const std::vector<mesh::FreezeMark>& marks);

    // The chord between two consecutive kept points of a border loop.
    struct KeptChord {
        float length = 0.0f;            // distance between the two points
        std::span<const int> loopEdges; // border edges of the loop; edge i leaves the point at position i
        int from = 0;                   // positions of the two points in the loop
        int to = 0;

        // Calls use(e) for the border edges the chord spans, in walking order.
        template <class Use>
        void forEachEdge(Use use) const
        {
            const int n = static_cast<int>(loopEdges.size());
            for (int i = from; i != to; i = (i + 1) % n) {
                use(loopEdges[static_cast<std::size_t>(i)]);
            }
        }
    };

    // Calls visit(chord) for the chords between consecutive kept points of
    // every border loop of the mesh: loop by loop, in each loop from one kept
    // position to the next in the order of keptPositions, the last chord back
    // to the first position.
    template <class Visit>
    static void forEachKeptChord(const mesh::Mesh& mesh, const std::vector<mesh::FreezeMark>& marks, Visit visit)
    {
        const std::vector<mesh::Vec3f>& points = mesh.points();
        const mesh::BoundaryLoops loops(mesh);
        for (int r = 0; r < loops.count(); ++r) {
            const auto loop = loops.points().row(r);
            const std::vector<int> kept = keptPositions(loop, marks);
            for (std::size_t k = 0; k < kept.size(); ++k) {
                KeptChord chord;
                chord.loopEdges = loops.edges().row(r);
                chord.from = kept[k];
                chord.to = kept[(k + 1) % kept.size()];
                chord.length = mesh::distance(points[static_cast<std::size_t>(loop[static_cast<std::size_t>(chord.to)])],
                                              points[static_cast<std::size_t>(loop[static_cast<std::size_t>(chord.from)])]);
                visit(chord);
            }
        }
    }
};

} // namespace areko::remesher::prep
