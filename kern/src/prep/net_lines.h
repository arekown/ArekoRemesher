#pragma once

#include "mesh/bit_array.h"
#include "mesh/csr.h"
#include "mesh/mesh.h"

#include <concepts>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace areko::remesher::prep {

// One line on the mesh: its points by number and position, a type and flags
// (1 from sharp edges, 2 face group border, 4 detected, 8 protection line).
struct NetLine {
    std::vector<int> vertices;
    std::vector<mesh::Vec3f> points;
    int type = 0;
    std::uint8_t flags = 0;
};

// The lines of the preparation and how they follow mesh changes.
//  - state -1: no sharp edge attribute, 0: the attribute is valid, 1: the list is valid,
//  - while the list is valid, the point map "line_origins" names for every
//    point its number when the list was last laid onto the mesh.
class NetLines {
public:
    static constexpr std::string_view originMap = "line_origins";
    // Point map (one bit per point) of the points on the final lines.
    static constexpr std::string_view pointMap = "line_points";

    enum class Restore { Direct = 1, OverTwoEdges = 2 };

    // Starts without lines: state 0 when the input brought sharp edges, else -1.
    void reset(bool hasEdgeAttribute);

    [[nodiscard]] int state() const noexcept { return state_; }
    [[nodiscard]] const std::vector<NetLine>& lines() const noexcept { return lines_; }
    // Hands the lines over and starts again without lines (state -1).
    [[nodiscard]] std::vector<NetLine> release();

    // Lines from the sharp edges (only with state 0 on a triangle mesh that
    // has sharp edges or face groups): point chains of the sharp edges with at
    // least two points become lines of type 2, flag 1; the list replaces the
    // old one, the origin map starts as identity, the sharp edges are removed
    // and the state becomes 1.
    void extractFromSharpEdges(mesh::Mesh& mesh);

    // Lays the lines onto the changed mesh over the origin map (the smallest
    // new point wins for an old one). A line starts at its first point that
    // still exists (none: it stays unchanged); following points that still
    // exist are appended, and their connection to the previous point is set
    // in `edges` when given: the direct edge, or with OverTwoEdges also a path
    // over a common neighbour, whose point is inserted into the line. Points
    // take the current positions; the origin map becomes identity again.
    void restore(mesh::Mesh& mesh, Restore mode, mesh::BitArray* edges);

    struct FinalLines {
        mesh::BitArray protectedEdges; // sharp edges and face group borders
        mesh::BitArray linePoints;     // one bit per point: on a line
    };

    // Edge sets the final lines are taken from.
    struct FinalSources {
        bool sharpEdges = true;   // the sharp edges (lines or attribute)
        bool groupBorders = true; // the face group borders
    };

    // Final lines, only on a triangle mesh and only when there is something
    // to take: with state 1 face groups or lines, with state 0 sharp edges,
    // otherwise face groups. First the chains of the sharp edges (with state 1
    // the edges of the lines restored onto the mesh, else the sharp edge
    // attribute) as type 2, flag 1, then the chains of the face group borders
    // (without the mesh border and without those sharp edges) as type 2,
    // flag 2, each only when `sources` names it; chains need two points. They
    // replace the list, the state becomes 1. Protected are both edge sets.
    [[nodiscard]] std::optional<FinalLines> buildFinal(mesh::Mesh& mesh, FinalSources sources);

    // Lines without flag bit 1 whose points `isDetected` names get flag bit 4.
    template <class Check>
        requires std::predicate<const Check&, std::span<const int>>
    void markDetected(const Check& isDetected)
    {
        for (NetLine& line : lines_) {
            if ((line.flags & 1u) == 0 && isDetected(std::span<const int>(line.vertices))) {
                line.flags = static_cast<std::uint8_t>(line.flags | 4u);
            }
        }
    }

    // One more line per chain (points in chain order at their current
    // positions) with the type of the chain and the given flags.
    void appendChains(const mesh::Mesh& mesh, const mesh::Csr& chains, std::span<const int> types,
                      std::uint8_t flags);

    // Lines with flags exactly 4 become protection lines (type 3, flag 8).
    void detectedToProtection();

    // Lines of type 2 with flag exactly 1 go back into the sharp edge
    // attribute (edge form, replaces it) and leave the list; state becomes 0.
    // Only with state 1.
    void toSharpEdges(mesh::Mesh& mesh);

    // Edge set of the preparation (one bit per edge of the current mesh):
    //  - groups: borders between face groups (with `withBoundary` also the
    //    mesh border), without face groups the border edges,
    //  - sharp: with state 1 the edges of the lines (restored onto the mesh
    //    with `restoreMode`), with state 0 the sharp edge attribute.
    [[nodiscard]] mesh::BitArray collectEdges(mesh::Mesh& mesh, bool groups, bool sharp, bool withBoundary,
                                              Restore restoreMode = Restore::Direct);

    // Lines that stand for sharp edges: type 2 with flag exactly 1.
    [[nodiscard]] static bool isSharpLine(const NetLine& line) noexcept { return line.type == 2 && line.flags == 1; }
    // Edges between consecutive points of the given lines (missing edges are skipped).
    static void addEdges(const mesh::Mesh& mesh, std::span<const NetLine> lines, bool onlySharp, mesh::BitArray& edges);

private:
    // How chain rows become lines.
    struct ChainForm {
        std::span<const int> types{}; // type per chain; empty: every line has type 2
        std::uint8_t flags = 0;
        bool validOnly = false;       // skip missing points (-1), drop lines with fewer than two points
    };

    void setIdentityOrigins(mesh::Mesh& mesh);
    // One line per chain row (points in chain order at their current
    // positions) in the given form, appended to `out`.
    static void appendChainLines(const mesh::Mesh& mesh, const mesh::Csr& chains, const ChainForm& form,
                                 std::vector<NetLine>& out);
    // Type 2 lines with the given flags from the point chains of an edge set
    // (chains need at least two points).
    [[nodiscard]] static std::vector<NetLine> linesFromEdges(const mesh::Mesh& mesh, const mesh::BitArray& edges,
                                                             std::uint8_t flags);

    int state_ = -1;
    std::vector<NetLine> lines_;
    int recordedPointCount_ = 0; // point count when the origin map was laid
};

} // namespace areko::remesher::prep
