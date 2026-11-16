#pragma once

#include "mesh/csr.h"
#include "mesh/mesh.h"

#include <vector>

namespace areko::remesher::extract {

// Points created during one round of the polygon splitting and their parents.
// Ids below the base count are points of the mesh, the others are new points
// in the order they were created. Every new point owns one row of parents
// (old point, weight); a parent that is itself new in this round is replaced
// by its own parents with the weights multiplied (float).
class NgonNewPoints {
public:
    explicit NgonNewPoints(const mesh::Mesh& mesh);

    [[nodiscard]] int base() const noexcept { return base_; }
    [[nodiscard]] int count() const noexcept { return static_cast<int>(positions_.size()); }
    [[nodiscard]] mesh::Vec3f position(int id) const;
    [[nodiscard]] const std::vector<mesh::Vec3f>& positions() const noexcept { return positions_; }
    [[nodiscard]] const mesh::Csr& parents() const noexcept { return parents_; }
    [[nodiscard]] const std::vector<float>& weights() const noexcept { return weights_; }

    // Appends a point and returns its id; its parents follow with addParent
    // and are closed with endParents.
    int add(mesh::Vec3f position);
    void addParent(int id, float weight);
    void endParents();
    // Moves a new point (id at or above base).
    void move(int id, mesh::Vec3f position);

    // Starts a new round on the mesh as it is now: base = point count, no new points.
    void reset();

private:
    const mesh::Mesh& mesh_;
    int base_;
    std::vector<mesh::Vec3f> positions_;
    mesh::Csr parents_;
    std::vector<float> weights_;
};

} // namespace areko::remesher::extract
