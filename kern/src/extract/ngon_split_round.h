#pragma once

#include "ngon_new_points.h"
#include "ngon_polygon_list.h"

#include "mesh/bit_array.h"
#include "mesh/csr.h"
#include "mesh/mesh.h"

#include <span>
#include <vector>

namespace areko::remesher::extract {

// One round of face replacements on a polygon mesh: faces marked as done are
// replaced by pending polygons, which may use new points. Commit applies the
// round to the mesh (kept faces in their order, then the pending polygons;
// point maps of new points follow their parents), reset starts the next round.
class SplitRound {
public:
    explicit SplitRound(mesh::Mesh& mesh);

    [[nodiscard]] mesh::Mesh& mesh() noexcept { return mesh_; }
    [[nodiscard]] NgonNewPoints& points() noexcept { return points_; }
    [[nodiscard]] int pendingCount() const noexcept { return pending_.rows(); }
    [[nodiscard]] bool done(int face) const noexcept { return done_.testBounded(static_cast<std::size_t>(face)); }

    // Splits face f (unless done) and queues its parts.
    void splitFace(int face);
    // Splits the polygon(s) of `list` standing for face f and queues the parts.
    void queueSplit(int face, NgonPolygonList& list);
    // Queues polygons that replace face f as they are.
    void queuePolygons(int face, const std::vector<std::vector<int>>& polygons);
    // All faces not done (the count of marks stays).
    void clearDone();

    void commit();
    void reset();

private:
    void markDone(int face);

    mesh::Mesh& mesh_;
    NgonNewPoints points_;
    mesh::Csr pending_;
    mesh::BitArray done_;
    int doneCount_ = 0;
};

} // namespace areko::remesher::extract
