#pragma once

#include "cross_field.h"
#include "curvature.h"

#include "mesh/bit_array.h"
#include "mesh/mesh.h"
#include "solver/sparse_matrix.h"

#include <array>
#include <vector>

namespace areko::remesher::field {

// The linear system of the cross field: one unknown angle per triangle that
// is computed, coupled across every edge whose quarter turn is settled, and
// held to the curvature direction at the anchors. Free edges and weak
// triangles enter later, when their integers are rounded.
class CrossFieldSystem {
public:
    // Weight of the coupling across an edge.
    static constexpr double smoothWeight = 200.0;
    // The quarter turn as float, widened: the unit of the integers in the system.
    static constexpr double quarterTurn = static_cast<double>(quarterTurnF);
    // Trust weight of the curvature from which a triangle is an anchor.
    static constexpr float anchorWeight = 500.0f;

    // What the system is built from; all of it outlives the system.
    struct Input {
        const mesh::Mesh& mesh;
        const Curvature& curvature;           // weights and the half that is not computed
        const mesh::BitArray& anchors;        // per triangle
        const mesh::BitArray& freeEdges;      // per edge
        const std::vector<float>& connection; // per edge
        const std::vector<double>& targetAngle;
    };

    // Angle of the curvature direction in the frame of every triangle
    // (-pi, pi]; triangles of the half that is not computed take pi minus the
    // angle of their mirror triangle (`faceMirror` null: none).
    [[nodiscard]] static std::vector<double> targetAngles(const mesh::Mesh& mesh, const Curvature& curvature,
                                                          const std::vector<mesh::Vec3f>& frameE1,
                                                          const std::vector<int>* faceMirror);

    // Quarter turns per edge that are settled before solving: across an inner
    // edge between two anchors (not both in the half that is not computed) the
    // nearest one between their curvature directions, 0 across every other edge.
    [[nodiscard]] static std::vector<int> anchorJumps(const Input& input);

    // Numbers the unknowns, lists the candidate edges and builds matrix and
    // right-hand side with the quarter turns `periodJump` across the edges
    // that are not free.
    CrossFieldSystem(const Input& input, const std::vector<int>& periodJump);

    [[nodiscard]] const Input& input() const noexcept { return input_; }
    // Unknown of a triangle, or -1 for a triangle that is not computed.
    [[nodiscard]] int unknownOf(int face) const { return unknownOfFace_[static_cast<std::size_t>(face)]; }
    [[nodiscard]] const std::vector<int>& unknownFaces() const noexcept { return unknownFaces_; }
    // Free edges between two computed triangles, ascending.
    [[nodiscard]] const std::vector<int>& candidateEdges() const noexcept { return candidates_; }
    // Unknowns of the two triangles of a candidate edge, in the order of the edge.
    [[nodiscard]] std::array<int, 2> edgeUnknowns(int edge) const;
    // A weight large enough to pull a triangle towards its curvature direction.
    [[nodiscard]] static bool carriesWeight(float weight) noexcept { return !(minimumWeight >= weight); }

    [[nodiscard]] solver::ListMatrix& matrix() noexcept { return matrix_; }
    [[nodiscard]] std::vector<double>& rhs() noexcept { return rhs_; }
    // Adds `value` to the entry of `column` in a sorted row, inserting it when absent.
    void add(int row, int column, double value);

private:
    static constexpr float minimumWeight = 2e-6f;

    void buildRows(const std::vector<int>& periodJump);

    Input input_;
    std::vector<int> unknownOfFace_;
    std::vector<int> unknownFaces_;
    std::vector<int> candidates_;
    solver::ListMatrix matrix_;
    std::vector<double> rhs_;
};

} // namespace areko::remesher::field
