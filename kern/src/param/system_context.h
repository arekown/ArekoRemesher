#pragma once

#include "parts/part.h"

#include <cstddef>
#include <vector>

namespace areko::remesher::param {

// What every step of the system set-up reads from a part: its mesh, the cut
// copy, the cut graph and the chains. The part has all of them and outlives
// the context.
class SystemContext {
public:
    // `symmetry`: the part may be computed by halves. Whether the part is
    // symmetric is asked once, here; that may look for its mirror pairs anew.
    SystemContext(parts::Part& part, bool symmetry)
        : part_(part), symmetry_(symmetry), pairs_(symmetry ? part.mirror.symmetric(part.mesh, true) : nullptr)
    {
    }

    [[nodiscard]] const parts::Part& part() const noexcept { return part_; }
    [[nodiscard]] const mesh::Mesh& mesh() const noexcept { return part_.mesh; }
    [[nodiscard]] const mesh::Adjacency& adjacency() const { return part_.mesh.adjacency(); }
    [[nodiscard]] const field::CutGraph& graph() const noexcept { return *part_.cutGraph; }
    [[nodiscard]] const field::CutCopy& copy() const noexcept { return part_.cutGraph->copy; }
    [[nodiscard]] const Chains& chains() const noexcept { return *part_.chains; }

    // Symmetry is switched on.
    [[nodiscard]] bool symmetry() const noexcept { return symmetry_; }
    // Symmetry is switched on and the part is symmetric.
    [[nodiscard]] bool symmetric() const noexcept { return pairs_ != nullptr; }
    // The mirror pairs of such a part, else null.
    [[nodiscard]] prep::SymmetryCache* pairs() const noexcept { return pairs_; }

    // Size of triangle `t` along the combed direction that belongs to grid
    // coordinate `coordinate` (0 or 1). A part without sizes per triangle has
    // its target edge length everywhere.
    [[nodiscard]] float sizeAlong(std::size_t t, int coordinate) const
    {
        if (!part_.sizeField) {
            return part_.edgeLength;
        }
        const field::SizeField& sizes = *part_.sizeField;
        const std::size_t faces = static_cast<std::size_t>(part_.mesh.faceCount());
        if (sizes.sizeU.size() != faces || sizes.sizeV.size() != faces) {
            return sizes.edgeLength;
        }
        return coordinate == 0 ? sizes.sizeU[t] : sizes.sizeV[t];
    }

    // The copy points of point `point` of the part mesh.
    [[nodiscard]] std::vector<int> copies(int point) const { return copy().pointCopies(adjacency(), point); }
    // The copy point of corner `corner`.
    [[nodiscard]] int copyPoint(int corner) const { return copy().topology().vertex(corner); }

    // Whether copy point `w` lies in the mirrored half: the triangle of its
    // first corner does.
    [[nodiscard]] bool mirrored(int w) const
    {
        const auto corners = copy().adjacency().pointCorners(w);
        return !corners.empty() && graph().mirrorHalf.testBounded(static_cast<std::size_t>(corners.front() / 3));
    }

private:
    const parts::Part& part_;
    bool symmetry_;
    prep::SymmetryCache* pairs_;
};

} // namespace areko::remesher::param
