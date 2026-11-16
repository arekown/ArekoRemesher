#pragma once

#include "trace.h"

#include "mesh/attribute_map.h"
#include "mesh/vec3.h"

#include <array>
#include <span>
#include <vector>

namespace areko::remesher::extract {

// Turns every closed walk of the trace into output polygons:
//  - a side whose segment came back to its own grid line is pulled together:
//    one of its two nodes is merged into the other (the end of a crease
//    curve, the start of the walk at a seed, the smaller number and the kept
//    points of a frozen border decide which one stays),
//  - a triangle left after merging without a node on a mesh point and with
//    a degenerate side still unmerged is dropped,
//  - a face of fewer than five corners (nodes taken over from a frozen
//    border not counted) becomes one polygon, unless it doubles back on
//    itself,
//  - a larger face first loses its spikes (a node that is met again two
//    steps on), then is cut into quads across its shortest diagonals (from a
//    corner to the third one after it), the rest becoming the last polygon.
class FaceClosing final : public FaceCloser {
public:
    // `frozenMarks`: the freeze marks of the part's points while its borders
    // are frozen, else null.
    explicit FaceClosing(const std::vector<mesh::FreezeMark>* frozenMarks) noexcept : frozenMarks_(frozenMarks) {}

    void close(FaceLoop& loop, Trace& trace) override;

    // Cuts the polygon `numbers` (corner positions `positions`, same order)
    // into quads: again and again the corner with the shortest distance to
    // the third corner after it (squared, float; the first one on a tie)
    // gives the quad of it and the next three corners, and the two corners
    // between are removed with their positions, until fewer than five are
    // left. The rest stays in `numbers`.
    [[nodiscard]] static std::vector<std::array<int, 4>> splitIntoQuads(std::vector<int>& numbers,
                                                                        std::vector<mesh::Vec3f> positions);

    // Output number of a node: the node it was merged into, else its own.
    [[nodiscard]] static int number(const TraceNode& node) noexcept;

private:
    // Merges the nodes of the degenerate sides; returns false when the face
    // is to be left out.
    bool merge(FaceLoop& loop, Trace& trace) const;
    static void emitSmall(const std::vector<int>& nodes, Trace& trace);
    // The nodes from `from` to `to` (taken modulo the node count, both
    // included) as one polygon; the indices between them go to `inner`.
    static void emitRange(const std::vector<int>& nodes, int from, int to, std::vector<bool>* inner, Trace& trace);
    static void emitSplit(const std::vector<int>& nodes, Trace& trace);
    static void endPolygon(Trace& trace);

    const std::vector<mesh::FreezeMark>* frozenMarks_;
};

} // namespace areko::remesher::extract
