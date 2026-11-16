#pragma once

#include "cost_queue.h"
#include "quadric.h"
#include "mesh/mesh.h"
#include "run/progress.h"

#include <cstdint>
#include <optional>
#include <vector>

namespace areko::remesher::prep {

// Edge-collapse decimation of a triangle mesh to a target number of vertices
// by plane quadrics (area weighted), single threaded with a fixed collapse
// order:
//  - vertex flags: free, line vertex (on a straight sharp edge or face group
//    border), locked (corner of such lines, three or more lines, border or
//    non-manifold); locked vertices never move, line vertices only along
//    their line,
//  - every vertex keeps its cheapest collapse; the globally cheapest one
//    (cost, vertex number) runs next after a recomputation and a fold-over
//    check; a failed collapse locks the vertex,
//  - colors and density factors are mixed by distance, face groups and sharp
//    edges follow.
// The source mesh is emptied; the result keeps only point colors, point
// density factors, face groups and sharp edges (per half edge).
class Decimator {
public:
    Decimator(float lengthScale, int targetVertexCount, bool useSharpEdges, run::Progress& progress);

    [[nodiscard]] mesh::Mesh run(mesh::Mesh&& source);

private:
    struct Evaluation {
        float cost = 0.0f;
        mesh::Vec3f position;
    };
    // Changes of one collapse, kept to undo it.
    struct Record {
        int vertex = -1;
        mesh::Vec3f oldPosition;
        std::vector<int> corners;
        std::vector<int> removedTriangles;
    };
    // New best collapse of a vertex, applied at the end of the round.
    struct Update {
        int vertex = -1;
        int target = -1;
        float cost = 0.0f;
        int bin = -1;
    };

    void load(mesh::Mesh& source);
    void computeFlags();
    void startCosts();
    void collapseLoop();
    [[nodiscard]] mesh::Mesh output() const;

    [[nodiscard]] static int nextCorner(int c) noexcept { return c + (c % 3 == 2 ? -2 : 1); }
    [[nodiscard]] static int prevCorner(int c) noexcept { return c + (c % 3 == 0 ? 2 : -1); }
    [[nodiscard]] std::vector<int> cornersOf(int v) const;
    [[nodiscard]] std::vector<int> neighbours(int v) const;
    [[nodiscard]] bool sharpBit(int halfEdge) const noexcept;
    void unlink(int vertex, int corner);
    void pushFront(int vertex, int corner);

    [[nodiscard]] Quadric triangleQuadric(int triangle) const;
    [[nodiscard]] Quadric vertexQuadric(const std::vector<int>& corners) const;
    [[nodiscard]] Evaluation evaluate(int a, int b, const Quadric& q) const;
    void fallback(int a, int b, const Quadric& q, bool lineCase, bool chooseB, float detAbs, Evaluation& e) const;
    [[nodiscard]] int check(int v, int t, mesh::Vec3f position, const std::vector<int>& cornersV,
                            const std::vector<int>& cornersT) const;
    [[nodiscard]] int collapse(int v, int t, mesh::Vec3f position, std::vector<int> cornersV, Record& record,
                               std::vector<int>& ringCorners, std::vector<int>& ring);
    void undo(int t, Record& record, std::uint8_t oldTargetFlag);
    void refresh(int t, const std::vector<int>& ringCorners, const std::vector<int>& ring, std::vector<Update>& updates);
    // Share of vertex t in a value mixed for `position` from the values of v
    // and t: the nearer vertex counts more.
    [[nodiscard]] float mixShare(int v, int t, mesh::Vec3f position) const;
    [[nodiscard]] mesh::Color4ub mixColor(int v, int t, float lambda) const;
    [[nodiscard]] float mixDensity(int v, int t, float lambda) const;

    // settings of the run
    float lengthLimit_ = 0.0f;   // cost below which long edges get a penalty
    float shortEdge2_ = 0.0f;    // squared length from which the penalty applies
    float penaltyOffset_ = 0.0f;
    float penaltyScale_ = 0.0f;
    int target_ = 0;
    bool useSharpEdges_ = false;
    run::Progress& progress_;

    // mesh state
    std::vector<int> tri_;            // vertex per corner
    std::vector<mesh::Vec3f> pos_;
    std::vector<mesh::Color4ub> colors_;
    bool mixColors_ = false;
    std::vector<float> density_;         // density factor per vertex, empty without
    std::vector<int> groups_;
    bool hasGroups_ = false;
    bool hasSharp_ = false;
    std::vector<std::uint8_t> triFlags_; // bits 0..2 sharp edge k -> k+1, bit 3 active
    std::vector<std::uint8_t> flags_;    // 1 active, 2 border, 4 locked, 8 line
    std::vector<int> head_;              // first corner per vertex
    std::vector<int> next_;              // next corner in the vertex list
    int active_ = 0;
    float diagonal_ = 0.0f;
    float weightScale_ = 0.0f;
    std::vector<Collapse> best_;
    std::vector<int> localOf_;          // scratch of refresh: local number per vertex, -1 = none
    std::vector<char> triangleMark_;    // scratch of refresh
    std::optional<CostQueue> queue_;
};

} // namespace areko::remesher::prep
