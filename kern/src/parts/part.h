#pragma once

#include "curves/guide_curve.h"
#include "extract/raw_mesh.h"
#include "extract/trace.h"
#include "field/cross_field.h"
#include "field/curvature.h"
#include "field/cut_graph.h"
#include "field/size_field.h"
#include "mesh/bit_array.h"
#include "mesh/mesh.h"
#include "param/chains.h"
#include "param/solution.h"
#include "param/system.h"
#include "prep/part_symmetry.h"
#include "prep/symmetry_cache.h"

#include <optional>
#include <vector>

namespace areko::remesher::parts {

// Mirror pairs of one part mesh and whether they are current. A real change
// of the mesh makes them old; the next query that may refresh looks for the
// symmetry anew (keeping a locked plane). Without a find the part is not
// symmetric; its point pairs and sides stay only while the point count is the
// same, all other pairs are dropped.
class PartMirror {
public:
    PartMirror() = default;
    // Pairs copied from the cut of the preparation (current, symmetric).
    explicit PartMirror(prep::SymmetryCache cache) : cache_(std::move(cache)), current_(true), symmetric_(true) {}

    // The mesh changed: the pairs are no longer current.
    void invalidate() noexcept { current_ = false; }

    // The pairs when the part is symmetric, else null. With `refresh`, old
    // pairs are looked for anew on `mesh` first.
    [[nodiscard]] prep::SymmetryCache* symmetric(const mesh::Mesh& mesh, bool refresh);
    // The stored pairs as they are (also when not symmetric), or null.
    [[nodiscard]] const prep::SymmetryCache* stored() const noexcept { return cache_ ? &*cache_ : nullptr; }

    // Evaluates the part against its pairs for a field step and stores the
    // result in `classes`. With `use`, old pairs are looked for anew on `mesh`
    // first and a symmetric part is evaluated; a part that is not symmetric,
    // and any part without `use` (symmetry switched off), gets the empty
    // evaluation. Returns the pairs of a symmetric part, else null.
    [[nodiscard]] prep::SymmetryCache* evaluate(const mesh::Mesh& mesh, bool use);

    // Evaluation of the part against its pairs, stored for the field steps.
    prep::PartSymmetry::Classes classes;

private:
    std::optional<prep::SymmetryCache> cache_;
    bool current_ = false;
    bool symmetric_ = false;
};

// One connected piece of the prepared mesh, computed on its own.
struct Part {
    int index = 0;
    mesh::Mesh mesh;                       // triangle mesh of the part
    std::vector<curves::GuideCurve> curves;
    int quadBudget = 0;
    float targetEdgeLength = -1.0f;        // given edge length (<= 0: none)
    float edgeLength = 0.0f;               // h: target edge length of the field steps
    bool skipped = false;                  // too small: no output
    PartMirror mirror;
    std::optional<field::Curvature> curvature; // empty when not computed or without radius
    std::optional<field::CrossField> crossField; // empty when not computed or not solvable
    std::optional<field::CutGraph> cutGraph;     // empty when not computed
    // Border edges whose length is not to shape the quads (one bit per edge):
    // the edges on the plane of a kept half. Empty when there is none.
    mesh::BitArray lockedBorderEdges;
    std::optional<field::SizeField> sizeField;   // empty when not computed
    std::optional<param::Chains> chains;         // empty when not computed
    std::optional<param::System> system;         // empty when not computed or not solvable
    std::optional<param::Solution> solution;     // empty when not computed or not solved
    std::optional<extract::Trace> trace;         // empty when not computed or nothing to trace
    std::optional<extract::RawMesh> raw;         // empty when not computed or nothing traced
    // Frozen borders give the part its sizes: it has freeze marks and the
    // mean chord between its kept border points is above 0. Set with the
    // target edge length. A part without border edge, or without a chord
    // that counts, has no border to freeze: it is computed like a part of a
    // mesh whose borders are not frozen, in every field step.
    bool bordersFrozen = false;

    // The freeze marks of the part's points while its borders are frozen, else null.
    [[nodiscard]] const std::vector<mesh::FreezeMark>* frozenBorderMarks() const;
};

} // namespace areko::remesher::parts
