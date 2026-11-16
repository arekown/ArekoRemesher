// Checks of core logic that can break silently: numbering rules of the
// neighbourhood, parts, settings conversions, map mixing, mesh changes, the
// order of equal keys in sorting, the bins of the cost queue, the regions of
// the group union, the mesh cleaning, edge flips (also required ones), the edge
// fan around a point, the crease limit in degrees, chains kept whole through
// break points, laying lines onto a changed mesh, the final lines, the repair runs of the
// smoother, the cut at a mirror plane, the detection of the symmetry, the
// split into parts with their quad budgets, the cut of a curve, the sparse
// solver (ordering, factorization and solving bit for bit, fixing, local update),
// the eigen decomposition of symmetric 3x3 matrices, the curvature of a cylinder,
// the singular points of the cross field and the cut into a disc on a sphere
// and a torus, the chains of the cuts with their nodes, sides and turns on
// both and on an open mesh, the start rules of the chains, a curve edge that
// gets no chain and a closed run of curve edges whose last piece joins its
// first.
// Also the elimination of linear conditions: column order, pivot rules, back
// substitution, arithmetic on term lists and grid periods, with values that
// are not numbers, infinite or out of range.
// Also the system of the parameterisation: variables, anchors, transitions
// for every quarter turn (every condition with its pivot and its terms in
// order, the rows before the elimination, rows that name a variable twice),
// an elimination that fails and leaves its rows to the system, a row term
// without an unknown, iso lines, energy rows and numbering on an open
// disc, with a cut, and on the sphere and the torus.
// Also its solution: the expected grid on a flat mesh, whole coordinates at
// nodes and transitions across the cuts, the order of the rounding steps,
// the reweighting with its averaging over neighbours, a stop request in a
// round, the rounding of what is left, the values after an elimination that
// failed, the mirror corners, the coordinates of a mirrored half and the
// start points.
// Also the rounding along frozen borders: what it selects and with which
// values, and where it stands in the order of the rounding steps, also
// behind the rounding by regions of the core; a point in two border loops,
// a condition term without unknown, a copy of the mirrored half, and the
// reference of a border curve, free or fixed before.
// Also the rounding by regions: the rows of the bridges for every kind of
// pair with their sides, signs and weights, two curves that meet, the size
// in the field sum, grid periods, and a run in which the regions find
// nothing left to deliver. And an elimination at the chain ends that fails.
// Also the repair of the raw mesh: every cleanup rule with its origin lists
// and with the maps kept, doublets at the border and inside, dissolving with
// a straight point and with a ring, split fans with their point data, holes,
// a joined thin hole (along both diagonals), notches and the gap rings of
// frozen borders.
// Also the trace of the whole-numbered grid lines over the triangles: nodes
// by kind and place, and every face of the grid walked once, on a flat
// square, on a plate around a point with three quads and across the seam
// of a cylinder; the fan the other way round around a shared point, the
// length limits of a walk, a second pass from the border points, and the
// place of a node in a triangle without numbers at one corner.
// Also the raw mesh: closing a face (spikes, merged sides, a dropped
// triangle), cutting a polygon into quads, the border classification by hand,
// and the whole raw mesh on the square, the plate around a point, the
// cylinder and the square with frozen borders: polygons by corner count,
// manifold, point sources after snapping onto corners.
// Also the result: giving the host groups back, splitting points along sharp
// edges, the maps and new face groups when parts come together, and the
// mirror at the end.
// Also the run: the parts on several threads give the result of one thread.
// Exit code 0 = all checks passed.

#include "areko_remesher.h"
#include "api/status.h"
#include "run/progress.h"
#include "solver/constrained_driver.h"
#include "solver/direct_solver.h"
#include "curves/curve_projector.h"
#include "curves/segment_cutter.h"
#include "curves/segment_preparation.h"
#include "field/border_lock.h"
#include "field/cross_field_solver.h"
#include "field/curvature_field.h"
#include "field/curve_alignment.h"
#include "field/cut_graph_builder.h"
#include "field/part_edge_length.h"
#include "field/size_field_builder.h"
#include "field/symmetric_eigen3.h"
#include "parts/part_builder.h"
#include "input/input_builder.h"
#include "mesh/attribute_maps.h"
#include "mesh/edge_chains.h"
#include "mesh/mesh.h"
#include "mesh/modification.h"
#include "mesh/sharp_edges.h"
#include "mesh/vertex_fan.h"
#include "param/anchor_constraints.h"
#include "param/chain_builder.h"
#include "param/chain_sides.h"
#include "param/constraint_periods.h"
#include "param/constraint_terms.h"
#include "param/curve_chain_builder.h"
#include "param/elimination.h"
#include "param/elimination_columns.h"
#include "param/transition_constraints.h"
#include "param/variable_flags.h"
#include "prep/cost_queue.h"
#include "prep/feature_chains.h"
#include "prep/flat_triangle_check.h"
#include "prep/group_union.h"
#include "prep/mesh_cleaner.h"
#include "prep/mesh_repair.h"
#include "prep/net_lines.h"
#include "prep/part_symmetry.h"
#include "prep/smoother.h"
#include "prep/symmetry_cut.h"
#include "prep/symmetry_detector.h"
#include "settings/settings.h"
#include "util/introsort.h"
#include "param/system_builder.h"
#include "param/reweighting.h"
#include "param/solution_builder.h"
#include "param/solution_seeds.h"
#include "param/uv_output.h"
#include "param/border_rounding.h"
#include "mesh/boundary_loops.h"
#include "param/bridge_rows.h"
#include "param/region_bridges.h"
#include "param/region_growth.h"
#include "param/region_rounding.h"
#include "param/seam_rows.h"
#include "param/soft_system.h"
#include "field/smoothing_neighbours.h"
#include "param/final_rounding.h"
#include "param/solve_strategy.h"
#include "mesh/degenerate_cleanup.h"
#include "extract/doublets.h"
#include "extract/edge_dissolve.h"
#include "extract/face_walker.h"
#include "extract/hole_closing.h"
#include "extract/isoline_tracer.h"
#include "extract/segment_tracer.h"
#include "extract/ngon_corner_classifier.h"
#include "extract/ngon_finish.h"
#include "extract/ngon_optimizer.h"
#include "extract/ngon_parity_balance.h"
#include "extract/ngon_polygon_splitter.h"
#include "extract/ngon_split_round.h"
#include "extract/ngon_strip_grid.h"
#include "extract/border_classifier.h"
#include "extract/extraction.h"
#include "extract/group_transfer.h"
#include "extract/hard_edge_transfer.h"
#include "extract/face_closing.h"
#include "extract/quad_mesh_builder.h"
#include "extract/final_cleanup.h"
#include "extract/final_relaxation.h"
#include "extract/frozen_border_reset.h"
#include "extract/source_weights.h"
#include "extract/surface_walker.h"
#include "extract/vertex_color_transfer.h"
#include "merge/crease_split.h"
#include "merge/group_restore.h"
#include "merge/result_assembly.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <functional>
#include <initializer_list>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using namespace areko::remesher;
using namespace areko::remesher::mesh;

namespace {

int failures = 0;
int checks = 0;

template <class A, class B>
void expect(const A& actual, const B& expected, const char* what, int line)
{
    ++checks;
    if (!(actual == expected)) {
        ++failures;
        std::printf("FEHLER Zeile %d: %s\n", line, what);
    }
}
#define EXPECT(actual, expected) expect((actual), (expected), #actual, __LINE__)

std::vector<int> v(std::span<const int> s) { return {s.begin(), s.end()}; }

Mesh polygonMesh(std::vector<Vec3f> points, std::vector<int> faceStart, std::vector<int> corners, FaceMode mode)
{
    Mesh m;
    m.swapPoints(points);
    m.editTopology().swapPolygons(faceStart, corners);
    if (mode != FaceMode::Polygons) {
        m.editTopology().trySetUniformMode(mode);
    }
    m.topologyChanged();
    return m;
}

std::vector<Vec3f> grid(int n)
{
    std::vector<Vec3f> p;
    for (int i = 0; i < n; ++i) {
        p.push_back({static_cast<float>(i), 0.0f, 0.0f});
    }
    return p;
}

void numberingTriangles()
{
    // t0 = (0,1,2), t1 = (2,1,3), t2 = (3,3,4) with a degenerate half edge 3 -> 3.
    Mesh m = polygonMesh(grid(5), {0, 3, 6, 9}, {0, 1, 2, 2, 1, 3, 3, 3, 4}, FaceMode::Triangles);
    EXPECT(m.topology().mode(), FaceMode::Triangles);
    const Adjacency& a = m.adjacency();
    EXPECT(a.edgeCount(), 6);
    EXPECT(a.cornerEdges(), (std::vector<int>{0, 2, 1, 2, 3, 4, -1, 5, 5}));
    EXPECT(a.edgeVertexList(), (std::vector<int>{0, 1, 0, 2, 1, 2, 1, 3, 2, 3, 3, 4}));
    EXPECT(v(a.edgeCorners(2)), (std::vector<int>{1, 3}));
    EXPECT(v(a.edgeCorners(5)), (std::vector<int>{7, 8}));
    EXPECT(v(a.pointCorners(3)), (std::vector<int>{5, 6, 7}));
    EXPECT(v(a.pointEdges(1)), (std::vector<int>{0, 2, 3}));
    EXPECT(v(a.pointEdges(3)), (std::vector<int>{3, 4, 5}));
    EXPECT(v(a.pointNeighbours(3)), (std::vector<int>{1, 2, 4}));
    EXPECT(a.boundaryEdges(), (std::vector<int>{0, 1, 3, 4}));
    EXPECT(a.boundaryPoints(), (std::vector<int>{0, 1, 2, 3}));
    EXPECT(a.next(2), 0);
    EXPECT(a.prev(0), 2);
    const ConnectedParts& parts = m.connectedParts();
    EXPECT(parts.partCount(), 2);
    EXPECT(parts.faceParts(), (std::vector<int>{0, 0, 1}));
    EXPECT(parts.partOfPoint(3), -2);
    EXPECT(parts.partOfPoint(4), 1);
}

void numberingQuads()
{
    // q0 = (0,1,4,3), q1 = (1,2,5,4).
    Mesh m = polygonMesh(grid(6), {0, 4, 8}, {0, 1, 4, 3, 1, 2, 5, 4}, FaceMode::Quads);
    EXPECT(m.topology().mode(), FaceMode::Quads);
    const Adjacency& a = m.adjacency();
    EXPECT(a.edgeCount(), 7);
    EXPECT(a.cornerEdges(), (std::vector<int>{0, 2, 5, 1, 3, 4, 6, 2}));
    EXPECT(a.edgeVertexList(), (std::vector<int>{0, 1, 0, 3, 1, 4, 1, 2, 2, 5, 3, 4, 4, 5}));
    EXPECT(v(a.edgeCorners(2)), (std::vector<int>{1, 7}));
    EXPECT(a.boundaryEdges(), (std::vector<int>{0, 1, 3, 4, 5, 6}));
    EXPECT(a.boundaryPoints(), (std::vector<int>{0, 1, 3, 2, 5, 4}));
    EXPECT(a.next(3), 0);
    EXPECT(a.prev(4), 7);
    EXPECT(a.counts().quads, 2);
    EXPECT(m.connectedParts().partCount(), 1);
}

void numberingPolygons()
{
    // f0 = pentagon (0,1,2,3,4), f1 = (1,0,5), f2 = (0,1,6): edge (0,1) has three corners.
    Mesh m = polygonMesh(grid(7), {0, 5, 8, 11}, {0, 1, 2, 3, 4, 1, 0, 5, 0, 1, 6}, FaceMode::Polygons);
    const Adjacency& a = m.adjacency();
    EXPECT(a.edgeCount(), 9);
    EXPECT(a.cornerEdges(), (std::vector<int>{0, 4, 7, 8, 1, 0, 2, 5, 0, 6, 3}));
    EXPECT(a.edgeVertexList(), (std::vector<int>{0, 1, 0, 4, 0, 5, 0, 6, 1, 2, 1, 5, 1, 6, 2, 3, 3, 4}));
    EXPECT(v(a.edgeCorners(0)), (std::vector<int>{0, 5, 8}));
    EXPECT(a.nonManifoldEdges(), (std::vector<int>{0}));
    EXPECT(a.boundaryEdges(), (std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8}));
    EXPECT(a.boundaryPoints(), (std::vector<int>{0, 4, 5, 6, 1, 2, 3}));
    EXPECT(v(a.pointEdges(0)), (std::vector<int>{0, 1, 2, 3}));
    EXPECT(a.edgeBetween(1, 0), 0);
    EXPECT(a.edgeBetween(4, 0), 1);
    EXPECT(a.edgeBetween(2, 4), -1);
    EXPECT(a.next(4), 0);
    EXPECT(a.prev(5), 7);
    EXPECT(a.counts().triangles, 2);
    EXPECT(a.counts().polygons, 1);
    EXPECT(m.connectedParts().partCount(), 1);
}

void conversions()
{
    using settings::SettingsConversion;
    float npc = 0, vdr = 0;
    const float cases[][3] = {{-10, 4, 1}, {0, 4, 1}, {25, 10, 2.5f}, {50, 16, 4}, {75, 24, 8}, {100, 32, 12}, {150, 32, 12}};
    for (const auto& c : cases) {
        SettingsConversion::curvatureToDensity(c[0], npc, vdr);
        EXPECT(npc, c[1]);
        EXPECT(vdr, c[2]);
    }
    SettingsConversion::curvatureToDensity(NAN, npc, vdr);
    EXPECT(std::isnan(npc) && std::isnan(vdr), true);
    EXPECT(SettingsConversion::quadCountFromPercentage(2048, 50.0f), 1024);
    EXPECT(SettingsConversion::quadCountFromPercentage(7872, 50.0f), 3936);
    EXPECT(SettingsConversion::quadCountFromPercentage(7872, 0.0f), 0);
    EXPECT(std::bit_cast<std::uint32_t>(SettingsConversion::degreesToRadians(30.0f)), 0x3f060a92u);
    EXPECT(SettingsConversion::inverseBorderFactor(2.0f), 0.5f);
    EXPECT(SettingsConversion::inverseBorderFactor(0.0f), 1.0f);

    arm_settings block;
    settings::Settings::fillDefaults(block);
    auto d = SettingsConversion::derive(block, 100);
    EXPECT(d.quadCount, 5000);
    EXPECT(d.pointsPerCircle, 16.0f);
    EXPECT(d.densityRatio, 4.0f);
    EXPECT(d.denoiseStrength, 0.5f);
    block.points_per_circle = 32.0f;
    block.quad_percentage = 50.0f;
    d = SettingsConversion::derive(block, 2048);
    EXPECT(d.pointsPerCircle, 32.0f);
    EXPECT(d.densityRatio, 4.0f);
    EXPECT(d.quadCount, 1024);
    block.quad_percentage = 0.0f;
    block.quad_count = 0;
    int status = ARM_OK;
    try {
        (void)SettingsConversion::derive(block, 10);
    } catch (const api::StatusError& e) {
        status = e.status();
    }
    EXPECT(status, ARM_ERROR_INVALID_SETTINGS);
    block.edge_length = 0.05f;
    EXPECT(SettingsConversion::derive(block, 10).edgeLength, 0.05f);

    settings::Settings::fillDefaults(block);
    block.mirror_x = 1;
    auto s = SettingsConversion::symmetry(block);
    EXPECT(s.on && s.mirrorAtEnd && s.curvatureSymmetry && s.axis[0] && !s.axis[1] && !s.halfMode, true);
    block.freeze_half = 1;
    s = SettingsConversion::symmetry(block);
    EXPECT(!s.on && !s.mirrorAtEnd && s.halfMode && s.axis[0], true);
    block.detect_symmetry = 1;
    s = SettingsConversion::symmetry(block);
    EXPECT(s.on && !s.mirrorAtEnd && !s.halfMode && !s.axis[0], true);
}

void mapMixing()
{
    Csr rows;
    rows.starts = {0, 2, 4, 4};
    rows.entries = {0, 1, 0, -1};
    TypedMap<int> ints("i", ElementKind::Point, {10, 20, -5});
    ints.interpolateAppend(rows, nullptr);
    EXPECT(ints.values(), (std::vector<int>{10, 20, -5, 15, 5, 0}));

    Csr one;
    one.starts = {0, 2};
    one.entries = {0, 1};
    const std::vector<float> w = {1.0f, 3.0f};
    TypedMap<int> weighted("i", ElementKind::Point, {10, 20});
    weighted.interpolateAppend(one, &w);
    EXPECT(weighted.values().back(), 17);

    TypedMap<std::uint8_t> bytes("b", ElementKind::Point, {10, 21});
    bytes.interpolateAppend(one, nullptr);
    EXPECT(bytes.values().back(), std::uint8_t{15});

    TypedMap<Color4ub> colors("c", ElementKind::Point, {Color4ub{1, 2, 3, 4}});
    colors.resize(2);
    EXPECT(colors.values()[1], (Color4ub{0, 0, 0, 0}));
    IndexMap renumber;
    renumber.newToOld = std::vector<int>{1, -1, 0};
    colors.reorder(renumber);
    EXPECT(colors.values()[1], (Color4ub{0, 0, 0, 255}));
    EXPECT(colors.values()[2], (Color4ub{1, 2, 3, 4}));

    BitArray bits(2);
    bits.set(0);
    BoolMap flags("f", ElementKind::Point, bits);
    Csr boolRows;
    // The last row names bit 2, which only this call appends (and sets): not a source.
    boolRows.starts = {0, 1, 3, 5, 6};
    boolRows.entries = {0, 0, 1, 0, -1, 2};
    flags.interpolateAppend(boolRows, nullptr);
    EXPECT(flags.bits().test(2), true);
    EXPECT(flags.bits().test(3), false);
    EXPECT(flags.bits().test(4), true);
    EXPECT(flags.bits().test(5), false);
}

void edgeChains()
{
    // Open chain 2-1-0 plus a separate closed loop 3-4-5: ends exist, branch points do not.
    Mesh m = polygonMesh(grid(6), {0, 3, 6}, {0, 1, 2, 3, 4, 5}, FaceMode::Triangles);
    const Adjacency& a = m.adjacency();
    BitArray chosen(static_cast<std::size_t>(a.edgeCount()));
    const std::pair<int, int> pairs[] = {{0, 1}, {1, 2}, {3, 4}, {4, 5}, {3, 5}};
    for (const auto& [p, q] : pairs) {
        chosen.set(static_cast<std::size_t>(a.edgeBetween(p, q)));
    }
    const Csr chains = EdgeChains::build(a, chosen);
    EXPECT(chains.starts, (std::vector<int>{0, 3, 7}));
    EXPECT(chains.entries, (std::vector<int>{2, 1, 0, 3, 4, 5, 3}));
}

// Strip of two quads: q0 = (0,1,4,3), q1 = (1,2,5,4).
Mesh strip()
{
    std::vector<Vec3f> p = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1.5f, 0}, {1, 1, 0}, {2, 1, 0}};
    return polygonMesh(p, {0, 4, 8}, {0, 1, 4, 3, 1, 2, 5, 4}, FaceMode::Polygons);
}

void setSharp(Mesh& m, std::initializer_list<std::pair<int, int>> pairs)
{
    BitArray bits(static_cast<std::size_t>(m.adjacency().edgeCount()));
    for (const auto& [a, b] : pairs) {
        bits.set(static_cast<std::size_t>(m.adjacency().edgeBetween(a, b)));
    }
    SharpEdges::setEdgeBits(m, std::move(bits));
}

bool isSharp(Mesh& m, int a, int b)
{
    const auto* edges = m.maps().findAs<BoolMap>(map_names::sharpEdges, ElementKind::Edge);
    const int e = m.adjacency().edgeBetween(a, b);
    return edges != nullptr && e >= 0 && static_cast<std::size_t>(e) < edges->bits().size() &&
           edges->bits().test(static_cast<std::size_t>(e));
}

void meshChanges()
{
    // Delete q0 and append the triangle (0,1,4) in one change.
    Mesh m = strip();
    m.maps().set(std::make_unique<TypedMap<int>>(std::string(map_names::materials), ElementKind::Face, std::vector<int>{7, 9}));
    m.maps().set(std::make_unique<TypedMap<float>>(std::string(map_names::density), ElementKind::Point,
                                                   std::vector<float>{1, 2, 3, 4, 5, 6}));
    setSharp(m, {{1, 4}, {0, 1}});
    Mesh triangle;
    {
        std::vector<int> fs = {0, 3}, cv = {0, 1, 4};
        triangle.editTopology().swapPolygons(fs, cv);
    }
    Modification change;
    change.deleteFaces = BitArray(2);
    change.deleteFaces->set(0);
    change.appendFaces = &triangle;
    change.points.state = 1; // the caller keeps the points
    EXPECT(MeshModifier::apply(m, change), true);
    EXPECT(m.topology().cornerVertex(), (std::vector<int>{1, 2, 5, 4, 0, 1, 4}));
    EXPECT(m.topology().faceStart(), (std::vector<int>{0, 4, 7}));
    EXPECT(m.maps().findAs<TypedMap<int>>(map_names::materials, ElementKind::Face)->values(), (std::vector<int>{9, 0}));
    EXPECT(m.maps().findAs<TypedMap<float>>(map_names::density, ElementKind::Point) != nullptr, true);
    EXPECT(isSharp(m, 1, 4), true);
    EXPECT(isSharp(m, 0, 1), false);
    EXPECT(change.halfEdges.newToOld->size(), std::size_t{7});

    // Pure deletion keeps point maps and sharp edges of the remaining faces.
    Mesh d = strip();
    setSharp(d, {{1, 4}});
    Modification removal;
    removal.deleteFaces = BitArray(2);
    removal.deleteFaces->set(1);
    EXPECT(MeshModifier::apply(d, removal), true);
    EXPECT(d.faceCount(), 1);
    EXPECT(removal.points.state, std::uint8_t{1});
    EXPECT(isSharp(d, 1, 4), true);

    // Complete new topology: sharp chains are laid onto the new faces; the
    // vanished diagonal (0,4) goes over the nearer common neighbour 1.
    std::vector<Vec3f> p = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1.5f, 0}, {1, 1, 0}, {2, 1, 0}};
    Mesh t = polygonMesh(p, {0, 3, 6, 10}, {0, 1, 4, 0, 4, 3, 1, 2, 5, 4}, FaceMode::Polygons);
    setSharp(t, {{0, 4}, {4, 5}});
    Modification replace;
    replace.faceStart = std::vector<int>{0, 4, 8};
    replace.cornerVertex = std::vector<int>{0, 1, 4, 3, 1, 2, 5, 4};
    replace.points.state = 1;
    replace.faces.state = 1; // same faces, new corner layout
    EXPECT(MeshModifier::apply(t, replace), true);
    EXPECT(isSharp(t, 0, 1) && isSharp(t, 1, 4) && isSharp(t, 4, 5), true);
    EXPECT(isSharp(t, 0, 3) || isSharp(t, 3, 4), false);
    EXPECT(*replace.faceStart, (std::vector<int>{0, 3, 6, 10}));
}

void replacementNormals()
{
    // t2 = (1,2,3) is degenerate (3 lies on 1-2) and takes the normal of t1.
    std::vector<Vec3f> p = {{0, 0, 0}, {2, 0, 0}, {0, 2, 0}, {1, 1, 0}};
    Mesh m = polygonMesh(p, {0, 3, 6, 9}, {0, 1, 3, 0, 3, 2, 1, 2, 3}, FaceMode::Triangles);
    EXPECT(m.triangleGeometry().normals()[2], (Vec3f{0, 0, 0}));
    EXPECT(m.triangleGeometry().totalArea(), 2.0f);
    m.setReplacementNormalMode(2, true);
    EXPECT(m.triangleGeometry().normals()[2], (Vec3f{0, 0, 1}));
    EXPECT(m.dihedralAngles()[static_cast<std::size_t>(m.adjacency().edgeBetween(0, 3))], 0.0f);
}

void inputBuilder()
{
    input::MeshInput in;
    const std::vector<float> points = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
    const std::vector<std::int32_t> sizes = {3, 3};
    const std::vector<std::int32_t> corners = {0, 1, 2, 0, 2, 3};
    in.setMesh(points, sizes, corners);
    // -0 equals +0; equal normals share one index.
    const std::vector<float> normals = {0, 0, 1, -0.0f, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, 1, 0, 0, 1};
    in.setCornerNormals(normals, true);
    const std::vector<float> colors = {1, 0, 0, 0.5f, 0, 1, 0, 1, 2.0f, -1.0f, 0, 1,
                                       0, 0, 1, 1, 0.5f, 0.5f, 0.5f, 1, 0, 0, 0, 1};
    in.setCornerColors(colors);
    arm_settings block;
    settings::Settings::fillDefaults(block);
    block.color_density = 1;
    Mesh m = input::InputBuilder::build(in, block);
    const auto& custom = *m.customNormals();
    EXPECT(custom.normals.size(), std::size_t{2});
    EXPECT(custom.cornerNormal[0], custom.cornerNormal[1]);
    EXPECT(custom.cornerNormal[3] != custom.cornerNormal[0], true);
    const auto* c = m.maps().findAs<TypedMap<Color4ub>>(map_names::colors, ElementKind::Point);
    EXPECT(c->values()[0], (Color4ub{0, 0, 255, 255}));   // last corner of point 0
    EXPECT(c->values()[2], (Color4ub{127, 127, 127, 255}));
    EXPECT(c->values()[1], (Color4ub{0, 255, 0, 255}));
}

// Order of equal keys after sorting (expected ids: the fixed course of steps
// of the introsort, including one case that ends in heap sort).
void introsortTies()
{
    using Item = std::pair<int, int>;
    auto sorted = [](const std::vector<int>& keys) {
        std::vector<Item> items;
        for (std::size_t i = 0; i < keys.size(); ++i) {
            items.push_back({keys[i], static_cast<int>(i)});
        }
        util::Introsort::sort(items.begin(), items.end(), [](const Item& a, const Item& b) { return a.first < b.first; });
        std::vector<int> ids;
        for (const Item& it : items) {
            ids.push_back(it.second);
        }
        return ids;
    };
    const std::vector<int> manyTies = {
        0, 0, 5, 5, 2, 1, 2, 6, 3, 5, 4, 4, 0, 3, 3, 0, 5, 2, 0, 1, 2, 0, 0, 1, 2, 2, 1, 3, 6, 3, 2, 6, 5, 1, 6, 6, 4,
        6, 5, 1, 3, 5, 4, 1, 5, 4, 1, 2, 2, 4, 6, 1, 0, 5, 3, 0, 4, 0, 3, 2, 5, 3, 6, 2, 5, 4, 2, 5, 3, 0, 1, 1, 3, 2,
        0, 6, 5, 1, 1, 6, 5, 1, 1, 6, 6, 3, 0, 4, 6, 3, 3, 3, 1, 3, 2, 4, 3, 4, 6, 0};
    EXPECT(sorted(manyTies),
           (std::vector<int>{1,  99, 86, 74, 69, 57, 55, 52, 22, 21, 18, 15, 12, 0,  33, 5,  51, 46, 43, 70,
                             39, 71, 77, 78, 81, 82, 26, 23, 19, 92, 6,  94, 59, 48, 47, 17, 63, 66, 4,  20,
                             73, 24, 25, 30, 72, 68, 93, 91, 13, 61, 85, 89, 58, 96, 90, 54, 14, 8,  27, 29,
                             40, 56, 97, 11, 95, 10, 87, 36, 42, 49, 45, 65, 41, 53, 60, 16, 9,  64, 2,  32,
                             67, 38, 76, 44, 3,  80, 75, 98, 7,  50, 37, 62, 28, 88, 35, 34, 84, 83, 31, 79}));
    // 17 elements: one partition step, then insertion.
    std::vector<int> seventeen;
    for (int i = 0; i < 17; ++i) {
        seventeen.push_back((i * 5) % 3);
    }
    EXPECT(sorted(seventeen), (std::vector<int>{0, 3, 15, 6, 9, 12, 2, 5, 14, 8, 11, 16, 13, 10, 7, 4, 1}));
    // Input that exhausts the depth limit (heap sort branch).
    const std::vector<int> deep = {23, 0,  26, 1,  30, 2,  25, 3,  28, 4,  24, 5,  31, 6,  23, 7,  27, 8,  22, 9,  29, 10,
                                   12, 11, 31, 30, 29, 28, 27, 26, 25, 24, 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                                   22, 20, 21, 19, 20, 18, 19, 17, 18, 16, 17, 15, 16, 14, 15, 13, 14, 12, 13, 21};
    EXPECT(sorted(deep), (std::vector<int>{1,  32, 3,  33, 5,  34, 7,  35, 9,  36, 11, 37, 13, 38, 15, 39,
                                           17, 40, 19, 41, 21, 42, 23, 43, 22, 61, 59, 62, 57, 60, 58, 55,
                                           56, 53, 51, 54, 49, 52, 47, 50, 45, 48, 63, 46, 18, 44, 14, 0,
                                           31, 10, 6,  30, 2,  29, 16, 28, 27, 8,  20, 26, 25, 4,  24, 12}));
}

// Bin number of the cost queue: logarithmic between the cost norm from the
// histogram and (diagonal / 100)^2, capped at the bin count; the front is the
// cheapest vertex by (cost, number).
void costQueueBins()
{
    const prep::FastLog2 log2;
    EXPECT(log2(1.0f), 0.0f);
    EXPECT(log2(8.0f), 3.0f);
    EXPECT(log2(0.5f), -1.0f);
    EXPECT(log2(1.5f), static_cast<float>(std::log(1.5) / 0.6931471824645996));

    // 100 vertices -> 10 bins; 20 costs of 1 set the norm to 10^0; diagonal 1000 -> upper end 10^2.
    std::vector<prep::Collapse> collapses(100, prep::Collapse{-1, std::numeric_limits<float>::max()});
    for (int v = 30; v < 50; ++v) {
        collapses[static_cast<std::size_t>(v)] = {0, 1.0f};
    }
    prep::CostQueue queue(collapses, 1000.0f);
    EXPECT(queue.binCount(), 10);
    EXPECT(queue.bin(0.5f), 0);
    EXPECT(queue.bin(1.0f), 0);
    EXPECT(queue.bin(3.0f), 2);
    EXPECT(queue.bin(50.0f), 8);
    EXPECT(queue.bin(1e6f), 10);
    EXPECT(queue.front(), 30);
    collapses[5] = {0, 50.0f};
    queue.insert(5, queue.bin(50.0f));
    for (int v = 30; v < 50; ++v) {
        queue.remove(v);
    }
    EXPECT(queue.queued(30), false);
    EXPECT(queue.front(), 5);
    collapses[7] = {0, 50.0f};
    queue.insert(7, queue.bin(50.0f));
    queue.remove(5);
    EXPECT(queue.front(), 7);
    queue.remove(7);
    EXPECT(queue.front(), -1);
}

// Group union with smoothing groups: regions are separated by material borders
// and smoothing borders; a smoothing border inside one region becomes sharp,
// one on a material border does not; existing sharp edges stay. Regions are
// numbered by their smallest face, the table holds the values of that face.
void groupUnionRegions()
{
    // 3 x 2 quads, faces q0 q1 q2 (bottom row) and q3 q4 q5 (top row).
    std::vector<Vec3f> p;
    for (int y = 0; y < 3; ++y) {
        for (int x = 0; x < 4; ++x) {
            p.push_back({static_cast<float>(x), static_cast<float>(y), 0.0f});
        }
    }
    Mesh m = polygonMesh(p, {0, 4, 8, 12, 16, 20, 24},
                         {0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 4, 5, 9, 8, 5, 6, 10, 9, 6, 7, 11, 10}, FaceMode::Quads);
    m.maps().set(std::make_unique<TypedMap<int>>(std::string(map_names::materials), ElementKind::Face,
                                                 std::vector<int>{0, 0, 1, 0, 0, 1}));
    // q0|q1 share no bit, q0-q3-q4-q1 do: one region around a smoothing border.
    m.maps().set(std::make_unique<TypedMap<std::uint32_t>>(std::string(map_names::smoothingGroups), ElementKind::Face,
                                                           std::vector<std::uint32_t>{1, 2, 4, 1, 3, 4}));
    setSharp(m, {{0, 1}});

    const prep::GroupRestoreTable table = prep::GroupUnion::apply(m);
    EXPECT(table.restore, true);
    EXPECT(table.materialsAreGroups, false);
    EXPECT(table.groupFaceGroup.empty(), true);
    EXPECT(table.groupMaterial, (std::vector<int>{0, 1}));
    EXPECT(table.groupSmoothing, (std::vector<std::uint32_t>{1, 4}));
    const auto* groups = m.maps().findAs<TypedMap<int>>(map_names::faceGroups, ElementKind::Face);
    EXPECT(groups != nullptr && groups->values() == (std::vector<int>{0, 0, 1, 0, 0, 1}), true);
    EXPECT(isSharp(m, 1, 5), true);
    EXPECT(isSharp(m, 0, 1), true);
    EXPECT(isSharp(m, 2, 6), false);
    EXPECT(isSharp(m, 6, 10), false);
    EXPECT(isSharp(m, 5, 6), false);
}

void meshCleaning()
{
    // Octahedron without its last face, plus a reversed copy of face 0: the
    // copy goes, the three point hole closes against the neighbours' direction.
    std::vector<Vec3f> p = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    Mesh m = polygonMesh(p, {0, 3, 6, 9, 12, 15, 18, 21, 24},
                         {0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 4, 2, 0}, FaceMode::Triangles);
    m.maps().set(std::make_unique<TypedMap<int>>(std::string(map_names::faceGroups), ElementKind::Face,
                                                 std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7}));
    prep::PrepState state;
    state.targetEdgeLength = 100.0f;
    prep::MeshCleaner::run(m, 0, state);
    EXPECT(m.topology().mode(), FaceMode::Triangles);
    EXPECT(m.faceCount(), 8);
    EXPECT(m.pointCount(), 6);
    EXPECT(m.adjacency().boundaryEdges().empty(), true);
    bool opposite = true;
    for (int e = 0; e < m.adjacency().edgeCount(); ++e) {
        const auto corners = m.adjacency().edgeCorners(e);
        opposite = opposite && corners.size() == 2 && m.adjacency().vertex(corners[0]) != m.adjacency().vertex(corners[1]);
    }
    EXPECT(opposite, true);
    EXPECT(m.maps().findAs<TypedMap<int>>(map_names::faceGroups, ElementKind::Face)->size(), std::size_t{8});
    EXPECT(state.repairCount, 0);
}

// The long diagonal of a thin rhombus is flipped; corners and the half edge
// form of the sharp edges follow, face values stay.
void flatTriangleFlip()
{
    std::vector<Vec3f> p = {{0, 0, 0}, {1, -0.2f, 0}, {2, 0, 0}, {1, 0.2f, 0}};
    Mesh m = polygonMesh(p, {0, 3, 6}, {0, 1, 2, 0, 2, 3}, FaceMode::Triangles);
    BitArray halfEdges(6);
    halfEdges.set(1); // 1 -> 2
    m.maps().set(std::make_unique<BoolMap>(std::string(map_names::sharpEdges), ElementKind::HalfEdge, halfEdges));
    m.maps().set(std::make_unique<TypedMap<int>>(std::string(map_names::faceGroups), ElementKind::Face, std::vector<int>{4, 7}));
    prep::FlatTriangleBlock block;
    prep::FlatTriangleCheck::run(m, BitArray(static_cast<std::size_t>(m.adjacency().edgeCount())), block);
    EXPECT(block.count, 1);
    EXPECT(block.again, true);
    EXPECT(m.topology().cornerVertex(), (std::vector<int>{0, 1, 3, 1, 2, 3}));
    const auto* sharp = m.maps().findAs<BoolMap>(map_names::sharpEdges, ElementKind::HalfEdge);
    EXPECT(sharp != nullptr && sharp->bits().count() == 1 && sharp->bits().test(3), true);
    const auto* groups = m.maps().findAs<TypedMap<int>>(map_names::faceGroups, ElementKind::Face);
    EXPECT((groups != nullptr && groups->values() == std::vector<int>{4, 7}), true);
    // Without a corner below 1 degree another run follows by round: mode 0
    // up to round 3, mode 3 in round 0, modes 1 and 2 never.
    const int modes[6] = {0, 0, 3, 3, 1, 2};
    const int rounds[6] = {3, 4, 0, 1, 0, 0};
    for (int i = 0; i < 6; ++i) {
        Mesh again = polygonMesh(p, {0, 3, 6}, {0, 1, 2, 0, 2, 3}, FaceMode::Triangles);
        prep::FlatTriangleBlock run;
        run.mode = modes[i];
        run.round = rounds[i];
        prep::FlatTriangleCheck::run(again, BitArray(static_cast<std::size_t>(again.adjacency().edgeCount())), run);
        EXPECT(run.count, 1);
        EXPECT(run.again, i == 0 || i == 2);
    }
}

void requiredFlatFlip()
{
    // t0 = (0,1,2) flat, t1 = (0,2,3) folded back by 135 degrees at 0-2; no
    // corner below 20 degrees. Mode 1 flips 0-2 only when 0 and 2 lie in the range.
    std::vector<Vec3f> p = {{0, 0, 0}, {1, -1, 0}, {2, 0, 0}, {1, -1, 1}};
    for (const int begin : {-1, 0}) {
        Mesh m = polygonMesh(p, {0, 3, 6}, {0, 1, 2, 0, 2, 3}, FaceMode::Triangles);
        prep::FlatTriangleBlock block;
        block.mode = 1;
        block.requiredBegin = begin;
        block.requiredEnd = 3;
        prep::FlatTriangleCheck::run(m, BitArray(static_cast<std::size_t>(m.adjacency().edgeCount())), block);
        EXPECT(block.count, begin < 0 ? 0 : 1);
        EXPECT(m.adjacency().edgeBetween(1, 3) >= 0, begin >= 0);
    }
}

void vertexFan()
{
    // Point 0 with the ring 1..4: t0 = (0,1,2), t1 = (0,2,3), t2 = (0,3,4), t3 = (0,4,1).
    std::vector<Vec3f> p = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0}, {0, -1, 0}};
    const Mesh closed = polygonMesh(p, {0, 3, 6, 9, 12}, {0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 1}, FaceMode::Triangles);
    const Adjacency& a = closed.adjacency();
    const int e1 = a.edgeBetween(0, 1), e2 = a.edgeBetween(0, 2), e3 = a.edgeBetween(0, 3), e4 = a.edgeBetween(0, 4);
    EXPECT(vertexFanEdges(a, 0, e1), (std::vector<int>{e1, e2, e3, e4, e1}));
    // Without t3 point 0 lies on the border: the fan runs from the outgoing border edge to the other one.
    const Mesh open = polygonMesh(p, {0, 3, 6, 9}, {0, 1, 2, 0, 2, 3, 0, 3, 4}, FaceMode::Triangles);
    const Adjacency& b = open.adjacency();
    EXPECT(vertexFanEdges(b, 0, b.edgeBetween(0, 1)),
           (std::vector<int>{b.edgeBetween(0, 1), b.edgeBetween(0, 2), b.edgeBetween(0, 3), b.edgeBetween(0, 4)}));
    EXPECT(vertexFanEdges(b, 0, b.edgeBetween(0, 4)).empty(), true);
}

// A ridge 0-1-2 folded by 80 degrees, point 1 raised so that the ridge bends
// by 3 degrees; each side is two flat strips (3-4-5 left, 6-7-8 right).
Mesh ridge()
{
    const float rise = std::tan(static_cast<float>(degD(1.5)));
    const float c = std::cos(static_cast<float>(degD(40)));
    const float s = std::sin(static_cast<float>(degD(40)));
    std::vector<Vec3f> p = {{0, 0, 0}, {1, 0, rise}, {2, 0, 0}};
    for (const float side : {c, -c}) {
        for (int i = 0; i < 3; ++i) {
            p.push_back({p[static_cast<std::size_t>(i)].x, side, p[static_cast<std::size_t>(i)].z - s});
        }
    }
    return polygonMesh(p, {0, 3, 6, 9, 12, 15, 18, 21, 24},
                       {0, 1, 4, 0, 4, 3, 1, 2, 5, 1, 5, 4, 1, 0, 6, 1, 6, 7, 2, 1, 7, 2, 7, 8}, FaceMode::Triangles);
}

void creaseLimitInDegrees()
{
    // A short chain with few neighbours: its smallest edge angle (80 degrees)
    // reaches the limit at L = 2 with h = 1.5 (75 degrees), it bends more than
    // its rivals: it stays a crease line.
    const Mesh m = ridge();
    prep::FeatureChains::Options options;
    options.lengths = {1.5f, 1.5f};
    const auto result = prep::FeatureChains::classify(m, options);
    EXPECT(result.classes, (std::vector<int>{2}));
}

void creaseChainBreak()
{
    // Point 0 is the corner of the flat quadrant 1-2-3 (z = 0); the rest of its
    // fan (ring 4 .. 11) falls away. The creases 0-1 (about 70 degrees) and 0-3
    // (about 63) turn by 90 degrees at 0 and jump by more than 5 degrees: 0 is a
    // break point. It does not split the chain 3-0-1, it only lowers its
    // smoothness, and the chain ends as a soft line.
    const float steep = std::tan(static_cast<float>(degD(60)));
    const float flat = std::tan(static_cast<float>(degD(50)));
    std::vector<Vec3f> p = {{0, 0, 0}, {-1, 0, 0}, {-1, -1, 0}, {0, -1, 0}};
    constexpr int ring = 8;
    for (int k = 1; k <= ring; ++k) {
        const double phi = degD(270.0 + 270.0 * k / (ring + 1));
        const float slope = flat + (steep - flat) * static_cast<float>(k) / static_cast<float>(ring + 1);
        p.push_back({static_cast<float>(std::cos(phi)), static_cast<float>(std::sin(phi)), -0.8f * slope});
    }
    std::vector<int> corners = {0, 1, 2, 0, 2, 3};
    std::vector<int> starts = {0, 3, 6};
    for (int k = 0, previous = 3; k <= ring; ++k) {
        const int next = k == ring ? 1 : 4 + k;
        corners.insert(corners.end(), {0, previous, next});
        starts.push_back(static_cast<int>(corners.size()));
        previous = next;
    }
    const Mesh m = polygonMesh(p, starts, corners, FaceMode::Triangles);
    prep::FeatureChains::Options options;
    options.creaseAngle = static_cast<float>(degD(50));
    options.lengths = {0.3f, 0.3f};
    const auto result = prep::FeatureChains::classify(m, options);
    EXPECT(result.points.entries, (std::vector<int>{3, 0, 1}));
    EXPECT(result.classes, (std::vector<int>{3}));
}

// 3 - 4 - 5     Lines from the sharp edges 0-1-2 are laid onto a mesh of the
// | / | \ |     same shape whose origin map names old 0 for the new points 3
// 0 - 1 - 2     and 5 (3 wins), old 1 for 4 and old 2 for 2. 4-2 is no edge:
// directly it is left out, over two edges it runs over 5 (shorter than over 1).
Mesh lineGrid()
{
    std::vector<Vec3f> p = {{0, 0, 0}, {0.5f, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}};
    return polygonMesh(p, {0, 3, 6, 9, 12}, {0, 1, 3, 1, 4, 3, 1, 2, 5, 1, 5, 4}, FaceMode::Triangles);
}

void lineRestore()
{
    for (const auto mode : {prep::NetLines::Restore::Direct, prep::NetLines::Restore::OverTwoEdges}) {
        Mesh before = lineGrid();
        setSharp(before, {{0, 1}, {1, 2}});
        prep::NetLines lines;
        lines.reset(true);
        lines.extractFromSharpEdges(before);
        EXPECT(lines.state(), 1);
        EXPECT(lines.lines().size(), std::size_t{1});
        if (lines.lines().size() != 1) {
            continue;
        }

        Mesh after = lineGrid();
        after.maps().set(std::make_unique<TypedMap<int>>(std::string(prep::NetLines::originMap), ElementKind::Point,
                                                         std::vector<int>{-1, -1, 2, 0, 1, 0}));
        BitArray edges;
        lines.restore(after, mode, &edges);
        const Adjacency& a = after.adjacency();
        const bool overTwo = mode == prep::NetLines::Restore::OverTwoEdges;
        std::vector<int> expected = overTwo ? std::vector<int>{3, 4, 5, 2} : std::vector<int>{3, 4, 2};
        std::vector<int> got = lines.lines().front().vertices;
        if (got.front() != expected.front()) {
            std::reverse(expected.begin(), expected.end());
        }
        EXPECT(got, expected);
        EXPECT(lines.lines().front().points.size(), expected.size());
        EXPECT(edges.size(), static_cast<std::size_t>(a.edgeCount()));
        EXPECT(edges.count(), std::size_t{overTwo ? 3u : 1u});
        EXPECT(edges.test(static_cast<std::size_t>(a.edgeBetween(3, 4))), true);
        if (overTwo) {
            EXPECT(edges.test(static_cast<std::size_t>(a.edgeBetween(4, 5))), true);
            EXPECT(edges.test(static_cast<std::size_t>(a.edgeBetween(5, 2))), true);
        }
    }
}

// 6 - 7 - 8    3 x 3 grid, two triangles per cell; the left faces form group 0,
// 3 - 4 - 5    the right ones group 1 (border 1-4-7). Sharp are 3-4 and 4-7:
// 0 - 1 - 2    the sharp line 3-4-7 comes first, the group line keeps only 1-4.
void finalLines()
{
    std::vector<Vec3f> p;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            p.push_back({static_cast<float>(c), static_cast<float>(r), 0.0f});
        }
    }
    std::vector<int> corners;
    std::vector<int> starts{0};
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 2; ++c) {
            const int a = 3 * r + c;
            corners.insert(corners.end(), {a, a + 1, a + 4, a, a + 4, a + 3});
            starts.push_back(static_cast<int>(corners.size()) - 3);
            starts.push_back(static_cast<int>(corners.size()));
        }
    }
    Mesh m = polygonMesh(p, starts, corners, FaceMode::Triangles);
    m.maps().set(std::make_unique<TypedMap<int>>(std::string(map_names::faceGroups), ElementKind::Face,
                                                 std::vector<int>{0, 0, 1, 1, 0, 0, 1, 1}));
    setSharp(m, {{3, 4}, {4, 7}});
    prep::NetLines lines;
    lines.reset(true);
    const auto result = lines.buildFinal(m, {});
    EXPECT(result.has_value(), true);
    EXPECT(lines.state(), 1);
    EXPECT(lines.lines().size(), std::size_t{2});
    if (!result || lines.lines().size() != 2) {
        return;
    }
    auto sorted = [](std::vector<int> x) {
        std::sort(x.begin(), x.end());
        return x;
    };
    EXPECT(lines.lines()[0].flags, std::uint8_t{1});
    EXPECT(sorted(lines.lines()[0].vertices), (std::vector<int>{3, 4, 7}));
    EXPECT(lines.lines()[1].flags, std::uint8_t{2});
    EXPECT(sorted(lines.lines()[1].vertices), (std::vector<int>{1, 4}));
    const Adjacency& a = m.adjacency();
    EXPECT(result->protectedEdges.count(), std::size_t{3});
    for (const auto& [u, w] : {std::pair{3, 4}, std::pair{4, 7}, std::pair{1, 4}}) {
        EXPECT(result->protectedEdges.test(static_cast<std::size_t>(a.edgeBetween(u, w))), true);
    }
    EXPECT(result->linePoints.size(), std::size_t{9});
    EXPECT(result->linePoints.count(), std::size_t{4});
    for (const int v : {1, 3, 4, 7}) {
        EXPECT(result->linePoints.test(static_cast<std::size_t>(v)), true);
    }
}

// Point 0 between t0 = (0,1,2) and t1 = (0,2,3) above the line 1-2-3; t0 is
// given the flipped normal, so the fold guard can satisfy neither triangle and
// puts point 0 back. After point 1 moved, the repair run finds a valid place.
void smootherRepair()
{
    std::vector<Vec3f> p = {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}, {1, 0, 0}};
    const Mesh m = polygonMesh(p, {0, 3, 6}, {0, 1, 2, 0, 2, 3}, FaceMode::Triangles);
    const std::vector<Vec3f> normals = {{0, 0, -1}, {0, 0, 1}};
    const std::vector<int> list = {1, 0};
    prep::Smoother::Options options;
    options.foldGuard = true;
    options.maxRepairRuns = -1;
    const std::vector<Vec3f> kept = prep::Smoother(m, options).applyTo(list, normals);
    EXPECT(kept[0], p[0]);
    options.maxRepairRuns = 4;
    const std::vector<Vec3f> repaired = prep::Smoother(m, options).applyTo(list, normals);
    EXPECT(repaired[1], kept[1]);
    EXPECT(repaired[0] == p[0], false);
    EXPECT(repaired[0].y > 0.1f && repaired[0].y < 0.25f && std::fabs(repaired[0].x) < 0.01f, true);
}

// Two triangles from x = -0.5 to 1.5 (plane through the middle, x = 0.5), the
// left border edge sharp, a point map with the old point numbers.
Mesh mirrorQuad()
{
    std::vector<Vec3f> p = {{-0.5f, 0, 0}, {1.5f, 0, 0}, {-0.5f, 1, 0}, {1.5f, 1, 0}};
    Mesh m = polygonMesh(p, {0, 3, 6}, {0, 1, 3, 0, 3, 2}, FaceMode::Triangles);
    setSharp(m, {{0, 2}});
    m.maps().set(std::make_unique<TypedMap<int>>("probe", ElementKind::Point, std::vector<int>{0, 1, 2, 3}));
    return m;
}

void mirrorCut()
{
    prep::SymmetryCut::Job job;
    job.splitEdgesFirst = true;
    job.usePresetPlane = true;
    Mesh m = mirrorQuad();
    const prep::SymmetryCut::Result result = prep::SymmetryCut::run(m, job);
    EXPECT(result.planePos, 0.5f);
    EXPECT(result.cache.has_value(), true);
    if (!result.cache) {
        return;
    }
    const auto& mirror = result.cache->pointMirror();
    const auto& points = m.points();
    EXPECT(mirror.size(), points.size());
    int onPlane = 0;
    bool mirrored = true;
    for (std::size_t i = 0; i < mirror.size(); ++i) {
        const auto j = static_cast<std::size_t>(mirror[i]);
        if (j == i) {
            ++onPlane;
            mirrored = mirrored && points[i].x == 0.5f;
            continue;
        }
        mirrored = mirrored && mirror[j] == static_cast<int>(i) && points[j].x == 1.0f - points[i].x &&
                   points[j].y == points[i].y && points[j].z == points[i].z &&
                   result.cache->pointSide().test(i) != result.cache->pointSide().test(j);
    }
    EXPECT(onPlane, 3);
    EXPECT(mirrored, true);
    const int faces = m.faceCount();
    const int h = faces / 2;
    EXPECT(faces, 6);
    bool facesMirrored = faces % 2 == 0;
    for (int f = 0; f < h; ++f) {
        facesMirrored = facesMirrored && result.cache->faceMirror()[static_cast<std::size_t>(f)] == h + f;
        for (int k = 0; k < 3; ++k) {
            const int a = m.topology().vertex(3 * f + k);
            const int b = m.topology().vertex(3 * (h + f) + (2 - k));
            facesMirrored = facesMirrored && mirror[static_cast<std::size_t>(a)] == b;
        }
    }
    EXPECT(facesMirrored, true);
    // Old points 0 and 2 are the first two kept ones.
    EXPECT(isSharp(m, 0, 1), true);
    EXPECT(isSharp(m, mirror[0], mirror[1]), true);

    // Keep-half: every kept point takes the map values of its own old point.
    Mesh half = mirrorQuad();
    job.cutOnly = true;
    (void)prep::SymmetryCut::run(half, job);
    const auto* probe = half.maps().findAs<TypedMap<int>>("probe", ElementKind::Point);
    EXPECT(probe != nullptr, true);
    if (probe != nullptr) {
        EXPECT(probe->values().size() >= 2 && probe->values()[0] == 0 && probe->values()[1] == 2, true);
    }
    EXPECT(half.faceCount(), 3);
}

// Point number of column x (-1 .. 1) in row y of mirrorStrip.
int stripPoint(int x, int y) { return 3 * y + x + 1; }

// A strip of 3 x 4 points mirrored at x = 0 (left cells cut by the diagonal
// from bottom left, right cells mirrored) in the plane z = 0, or z = y when
// `tilted`, with the middle column moved along x by `shift` per row. The
// points `extra` and the triangles `extraCorners` follow the strip.
Mesh mirrorStrip(std::array<float, 4> shift, bool tilted, std::initializer_list<Vec3f> extra = {},
                 std::initializer_list<int> extraCorners = {})
{
    std::vector<Vec3f> p;
    for (int y = 0; y < 4; ++y) {
        const float fy = static_cast<float>(y);
        for (int x = -1; x <= 1; ++x) {
            const float fx = static_cast<float>(x) + (x == 0 ? shift[static_cast<std::size_t>(y)] : 0.0f);
            p.push_back({fx, fy, tilted ? fy : 0.0f});
        }
    }
    p.insert(p.end(), extra);
    std::vector<int> corners;
    for (int y = 0; y < 3; ++y) {
        corners.insert(corners.end(),
                       {stripPoint(-1, y), stripPoint(0, y), stripPoint(0, y + 1), stripPoint(-1, y),
                        stripPoint(0, y + 1), stripPoint(-1, y + 1), stripPoint(1, y), stripPoint(0, y + 1),
                        stripPoint(0, y), stripPoint(1, y), stripPoint(1, y + 1), stripPoint(0, y + 1)});
    }
    corners.insert(corners.end(), extraCorners);
    std::vector<int> starts;
    for (std::size_t c = 0; c <= corners.size(); c += 3) {
        starts.push_back(static_cast<int>(c));
    }
    return polygonMesh(std::move(p), std::move(starts), std::move(corners), FaceMode::Triangles);
}

// mirrorStrip and two separate triangles above it at mirrored places, which
// only the distance pairs can match.
void mirrorDetection()
{
    const Mesh m = mirrorStrip({}, false,
                               {{-0.7f, 1.1f, 1.0f}, {-0.3f, 1.2f, 1.0f}, {-0.5f, 1.7f, 1.0f},
                                {0.7f, 1.1f, 1.0f}, {0.3f, 1.2f, 1.0f}, {0.5f, 1.7f, 1.0f}},
                               {12, 13, 14, 15, 17, 16});
    const std::vector<Vec3f>& p = m.points();
    std::optional<prep::SymmetryCache> cache = prep::SymmetryDetector::run(m, nullptr);
    EXPECT(cache.has_value(), true);
    if (!cache) {
        return;
    }
    const auto& mirror = cache->pointMirror();
    bool points = mirror.size() == p.size();
    for (std::size_t i = 0; points && i < mirror.size(); ++i) {
        const int j = mirror[i];
        points = j >= 0 && mirror[static_cast<std::size_t>(j)] == static_cast<int>(i) &&
                 p[static_cast<std::size_t>(j)].x == -p[i].x && p[static_cast<std::size_t>(j)].y == p[i].y &&
                 p[static_cast<std::size_t>(j)].z == p[i].z && (p[i].x != 0.0f || j == static_cast<int>(i));
        // Side 1 lies towards the normal (+x).
        points = points && (i >= 12 || p[i].x == 0.0f || cache->pointSide().test(i) == (p[i].x > 0.0f));
    }
    EXPECT(points, true);
    bool faces = true;
    for (int f = 0; f < 12; ++f) {
        const int g = cache->faceMirror()[static_cast<std::size_t>(f)];
        const bool right = (f % 4) >= 2;
        faces = faces && g >= 0 && g < 12 && cache->faceMirror()[static_cast<std::size_t>(g)] == f &&
                cache->faceSide().test(static_cast<std::size_t>(f)) == right &&
                cache->sideOfFace(m.topology(), f) == (right ? 1 : 0);
        for (int k = 0; faces && k < 3; ++k) {
            const int a = mirror[static_cast<std::size_t>(m.topology().vertex(3 * f + k))];
            faces = m.topology().vertex(3 * g) == a || m.topology().vertex(3 * g + 1) == a ||
                    m.topology().vertex(3 * g + 2) == a;
        }
    }
    EXPECT(faces, true);
    // The two separate triangles keep no face pair.
    EXPECT(cache->faceMirror()[12] == -1 && cache->faceMirror()[13] == -1, true);
    EXPECT(cache->mirroredFaceCount(), 12);
    EXPECT(cache->allFacesMirrored(m), false);
    const Adjacency& a = m.adjacency();
    EXPECT(cache->planeEdges(a),
           (std::vector<int>{a.edgeBetween(stripPoint(0, 0), stripPoint(0, 1)),
                             a.edgeBetween(stripPoint(0, 1), stripPoint(0, 2)),
                             a.edgeBetween(stripPoint(0, 2), stripPoint(0, 3))}));
    const prep::PartSymmetry::Classes classes = prep::PartSymmetry::classify(m, *cache);
    EXPECT(classes.valid && classes.mirrorGroup == 0, true);
    const BitArray half = prep::PartSymmetry::half(m, *cache, classes);
    bool left = half.size() == 14;
    for (int f = 0; left && f < 14; ++f) {
        left = half.test(static_cast<std::size_t>(f)) == (f < 12 && (f % 4) < 2);
    }
    EXPECT(left, true);
}

// mirrorStrip tilted into the plane z = y, with the middle column moved along
// x per row. A line that leaves the plane x = 0 must stay straight in the
// strip: a turn of about 4 degrees keeps the line and the symmetry, one of
// about 19 degrees drops the line, and with it the only start edge.
void mirrorLineBend()
{
    const Mesh slight = mirrorStrip({0.0f, 0.0f, 0.1f, 0.2f}, true);
    std::optional<prep::SymmetryCache> kept = prep::SymmetryDetector::run(slight, nullptr);
    EXPECT(kept.has_value() && kept->planeEdges(slight.adjacency()).size() == 3, true);
    const Mesh bent = mirrorStrip({0.0f, 0.0f, 0.5f, 0.5f}, true);
    EXPECT(prep::SymmetryDetector::run(bent, nullptr).has_value(), false);
}

} // namespace

// Squares of the given sides (two triangles each) side by side, apart.
Mesh squares(std::initializer_list<float> sides)
{
    std::vector<Vec3f> p;
    std::vector<int> corners;
    float x = 0.0f;
    for (const float s : sides) {
        const int b = static_cast<int>(p.size());
        p.insert(p.end(), {{x, 0.0f, 0.0f}, {x + s, 0.0f, 0.0f}, {x + s, s, 0.0f}, {x, s, 0.0f}});
        corners.insert(corners.end(), {b, b + 1, b + 2, b, b + 2, b + 3});
        x += s + 1.0f;
    }
    Mesh m;
    m.swapPoints(p);
    m.editTopology().swapCorners(corners, FaceMode::Triangles);
    m.topologyChanged();
    return m;
}

void partsAndBudget()
{
    // Areas 1 : 2 : 1000; a face map goes along.
    Mesh m = squares({1.0f, std::sqrt(2.0f), std::sqrt(1000.0f)});
    m.maps().set(std::make_unique<TypedMap<int>>(std::string(map_names::materials), ElementKind::Face,
                                                 std::vector<int>{0, 1, 2, 3, 4, 5}));
    parts::PartBuilder::Split split = parts::PartBuilder::split(m, nullptr, false);
    EXPECT(split.sharedPoints, false);
    EXPECT(split.parts.size(), std::size_t{3});
    if (split.parts.size() != 3) {
        return;
    }
    EXPECT(split.parts[2].mesh.points()[0] == m.points()[8], true);
    EXPECT(v(split.parts[1].mesh.topology().cornerVertex()), (std::vector<int>{0, 1, 2, 0, 2, 3}));
    EXPECT(split.facePlaces, (std::vector<int>{0, 1, 0, 1, 0, 1}));
    const auto* partMaterials = split.parts[2].mesh.maps().findAs<TypedMap<int>>(map_names::materials, ElementKind::Face);
    EXPECT(partMaterials != nullptr && partMaterials->values() == (std::vector<int>{4, 5}), true);
    parts::PartBuilder::assignBudgets(split.parts, m, 5000, -1.0f);
    // A_i * (N / A) in float, at least 10, truncated.
    const float total = m.triangleGeometry().totalArea();
    const float perArea = 5000.0f / total;
    EXPECT(split.parts[0].quadBudget, 10);
    EXPECT(split.parts[1].quadBudget, 10);
    EXPECT(split.parts[2].quadBudget, truncToInt(split.parts[2].mesh.triangleGeometry().totalArea() * perArea));
    EXPECT(split.parts[2].quadBudget, 4985);
    // One part: trunc(N) without the lower limit; with an edge length N = A / L^2.
    Mesh one = squares({1.0f});
    parts::PartBuilder::Split single = parts::PartBuilder::split(one, nullptr, false);
    parts::PartBuilder::assignBudgets(single.parts, one, 5000, 0.5f);
    EXPECT(single.parts.size() == 1 && single.parts[0].quadBudget == 4, true);
    parts::PartBuilder::assignBudgets(single.parts, one, 7, -1.0f);
    EXPECT(single.parts[0].quadBudget, 7);
    // A square (two triangles on four points) is not skipped.
    EXPECT(parts::PartBuilder::negligible(single.parts[0], m), false);
}

void curveCut()
{
    // Strip of 4 x 3 squares (two triangles each); a curve across it at y = 1.5.
    std::vector<Vec3f> p;
    for (int j = 0; j <= 3; ++j) {
        for (int i = 0; i <= 4; ++i) {
            p.push_back({static_cast<float>(i), static_cast<float>(j), 0.0f});
        }
    }
    std::vector<int> corners;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 4; ++i) {
            const int a = j * 5 + i;
            corners.insert(corners.end(), {a, a + 1, a + 6, a, a + 6, a + 5});
        }
    }
    Mesh m;
    m.swapPoints(p);
    m.editTopology().swapCorners(corners, FaceMode::Triangles);
    m.topologyChanged();
    curves::LocatedPolyline line;
    line.points = {{0.25f, 1.5f, 0.0f}, {3.75f, 1.5f, 0.0f}};
    line.locations = {{curves::Location::Face, 9}, {curves::Location::Face, 14}};
    const std::vector<curves::CurveSegment> pieces = curves::CurveProjector::project(m, line);
    EXPECT(pieces.size(), std::size_t{1});
    std::vector<const curves::CurveSegment*> list;
    for (const auto& piece : pieces) {
        list.push_back(&piece);
    }
    const std::vector<curves::CurveSegment> prepared =
        curves::SegmentPreparation::prepare(m, list, curves::SegmentPreparation::Options{});
    const curves::SegmentCutter::Result result = curves::SegmentCutter::cut(m, prepared, curves::SegmentCutter::Options{});
    EXPECT(result.changed, true);
    EXPECT(result.chains.size(), std::size_t{1});
    if (result.chains.size() != 1) {
        return;
    }
    const std::vector<int>& chain = result.chains[0];
    bool joined = chain.size() >= 2;
    bool onLine = true;
    for (std::size_t i = 0; i < chain.size(); ++i) {
        onLine = onLine && std::fabs(m.points()[static_cast<std::size_t>(chain[i])].y - 1.5f) < 1e-6f;
        if (i > 0) {
            joined = joined && m.adjacency().edgeBetween(chain[i - 1], chain[i]) >= 0;
        }
    }
    EXPECT(joined, true);
    EXPECT(onLine, true);
    EXPECT(m.faceCount() > 24, true);
    EXPECT(m.adjacency().nonManifoldEdges().empty(), true);
}

namespace {

// Small symmetric system in full storage (columns ascending) with its right-hand side.
struct SolverCase {
    int n;
    std::vector<int> starts;
    std::vector<int> rows;
    std::vector<double> values;
    std::vector<double> rhs;
};

// Expected values computed step by step (ordering, factorization with division
// by D, solving with the reciprocal of D).
const SolverCase gridCase{9,
    {0, 3, 7, 10, 14, 19, 23, 26, 30, 33},
    {0, 1, 3, 0, 1, 2, 4, 1, 2, 5, 0, 3, 4, 6, 1, 3, 4, 5, 7, 2, 4, 5, 8, 3, 6, 7, 4, 6, 7, 8, 5, 7, 8},
    {4.25, -1.0, -0.9, -1.0, 4.375, -1.0625, -0.87, -1.0625, 4.5, -0.8400000000000001, -0.9, 4.625, -1.1875,
     -0.81, -0.87, -1.1875, 4.75, -1.25, -0.78, -0.8400000000000001, -1.25, 4.875, -0.75, -0.81, 5.0, -1.375,
     -0.78, -1.375, 5.125, -1.4375, -0.75, -1.4375, 5.25},
    {0.3333333333333333, 0.8333333333333333, 1.7333333333333332, 2.233333333333333, 3.1333333333333333,
     3.6333333333333333, 4.533333333333332, 5.033333333333332, 5.933333333333333}};
const std::vector<int> gridListOrder{8, 7, 6, 3, 4, 0, 2, 5, 1};
const std::vector<std::uint64_t> gridListBits{
    0x3fe37cc1ecd5109bull, 0x3feeca27bd9e0248ull, 0x3fedfc595e2ee1c8ull, 0x3ff6faffa7079558ull,
    0x40004a5e75b9e17full, 0x3ffbd49783b83fa9ull, 0x3ffc823b7859f7a9ull, 0x4002b08242e9a4e3ull,
    0x40002582327c4567ull};
const SolverCase arrowCase{7,
    {0, 3, 7, 11, 17, 21, 24, 25},
    {0, 1, 3, 0, 1, 2, 3, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 2, 3, 4, 5, 3, 4, 5, 6},
    {3.1, 0.8, 0.35, 0.8, -2.5, -1.1, 0.44999999999999996, -1.1, 4.7, 0.55, 0.6, 0.35, 0.44999999999999996, 0.55,
     9.3, -0.25, -0.2, 0.6, -0.25, 3.1, -0.7, -0.2, -0.7, 3.3, 3.7},
    {0.3, -1.7, 2.2, 0.9, -0.4, 1.3, 0.1}};
const std::vector<int> arrowListOrder{5, 4, 2, 0, 1, 6, 3};
const std::vector<std::uint64_t> arrowListBits{
    0xbf92d7f09a62aeccull, 0x3fdb4650eda18e09ull, 0x3fe2a28e42b4fb96ull, 0x3fa791f3fb9897c0ull,
    0xbfc3f508ac7eba5dull, 0x3fd74625c6e242b8ull, 0x3f9bacf914c1bacfull};
const SolverCase zeroCase{8,
    {0, 3, 7, 11, 16, 21, 25, 29, 32},
    {0, 3, 4, 1, 3, 4, 5, 2, 3, 4, 6, 0, 1, 2, 3, 5, 0, 1, 2, 4, 7, 1, 3, 5, 6, 2, 5, 6, 7, 4, 6, 7},
    {5.5, -1.2, 0.7, 5.25, -0.8, 0.9, 0.0, 5.0, 1.1, -0.6, 0.0, -1.2, -0.8, 1.1, 4.75, 0.4, 0.7, 0.9, -0.6, 4.5,
     -0.35, 0.0, 0.4, 4.25, -1.3, 0.0, -1.3, 4.0, 0.55, -0.35, 0.55, 3.75},
    {1.1, -0.3, 0.75, 2.0, -1.4, 0.6, 0.2, -0.9}};
const std::vector<int> zeroListOrder{7, 6, 2, 0, 1, 4, 5, 3};
const std::vector<std::uint64_t> zeroListBits{
    0x3fd76e855534ad16ull, 0x3fb7c5a23815598eull, 0xbf8cd2eb6efcb1b8ull, 0x3fe0ade7ec007904ull,
    0xbfda594598e567e6ull, 0x3fc10c7a644e6296ull, 0x3fc1301007ed25fdull, 0xbfd31463d4d32ce9ull};
const std::vector<int> zeroDroppedOrder{7, 6, 5, 0, 1, 2, 4, 3};
const std::vector<std::uint64_t> zeroDroppedBits{
    0x3fd76e855534ad15ull, 0x3fb7c5a23815598dull, 0xbf8cd2eb6efcb1a8ull, 0x3fe0ade7ec007903ull,
    0xbfda594598e567e6ull, 0x3fc10c7a644e6294ull, 0x3fc1301007ed25fcull, 0xbfd31463d4d32ce9ull};
const std::vector<std::uint64_t> greedyFixBits{
    0x3fe0d0d0d0d0d0d1ull, 0x3ff0000000000000ull, 0x3fee13fdf1460886ull, 0x3ff0000000000000ull,
    0x3ffeb13108c37884ull, 0x3ffb4ff66cdf250dull, 0x3ffb2e0bf24566c8ull, 0x4002545781424bccull,
    0x400002cc81eca810ull};
const std::vector<std::uint64_t> cgUpdateBits{
    0xbf92d7f09a62aeccull, 0x3fd354cee9e0098bull, 0x3feb381f8150ef88ull, 0x3fa1e795fc13d6d5ull,
    0xbfd525f23678b79cull, 0x3fd6c5a07004c829ull, 0x3f9bacf914c1bacfull};

solver::CompactMatrix compactOf(const SolverCase& c)
{
    solver::CompactMatrix m(c.n);
    for (int col = 0; col < c.n; ++col) {
        for (int k = c.starts[static_cast<std::size_t>(col)]; k < c.starts[static_cast<std::size_t>(col) + 1]; ++k) {
            m.append(c.rows[static_cast<std::size_t>(k)], c.values[static_cast<std::size_t>(k)]);
        }
        m.closeColumn();
    }
    return m;
}

solver::ListMatrix listOf(const SolverCase& c)
{
    solver::ListMatrix m(c.n, c.n);
    for (int col = 0; col < c.n; ++col) {
        for (int k = c.starts[static_cast<std::size_t>(col)]; k < c.starts[static_cast<std::size_t>(col) + 1]; ++k) {
            m.list(col).push_back({c.rows[static_cast<std::size_t>(k)], c.values[static_cast<std::size_t>(k)]});
        }
    }
    return m;
}

std::vector<std::uint64_t> bitsOf(const std::vector<double>& x)
{
    std::vector<std::uint64_t> out;
    for (const double v : x) {
        out.push_back(std::bit_cast<std::uint64_t>(v));
    }
    return out;
}

struct QuietProgress {
    std::atomic<bool> cancel{false};
    run::Progress progress{nullptr, nullptr, cancel};
};

void directSolverOrders()
{
    QuietProgress q;
    solver::DirectSolver solver;
    std::vector<double> x;
    // Equal degrees on a grid, supernodes and mass elimination.
    EXPECT(solver.solve(listOf(gridCase), gridCase.rhs, q.progress, x), true);
    EXPECT(solver.ordering(), gridListOrder);
    EXPECT(bitsOf(x), gridListBits);
    EXPECT(solver.solve(compactOf(gridCase), false, gridCase.rhs, q.progress, x), true);
    EXPECT(bitsOf(x), gridListBits);
    // Dense node last, isolated node, negative D.
    EXPECT(solver.solve(compactOf(arrowCase), true, arrowCase.rhs, q.progress, x), true);
    EXPECT(solver.ordering(), arrowListOrder);
    EXPECT(bitsOf(x), arrowListBits);
    // Explicit zeros count in the list path and are dropped in the compressed path.
    EXPECT(solver.solve(listOf(zeroCase), zeroCase.rhs, q.progress, x), true);
    EXPECT(solver.ordering(), zeroListOrder);
    EXPECT(bitsOf(x), zeroListBits);
    EXPECT(solver.solve(compactOf(zeroCase), true, zeroCase.rhs, q.progress, x), true);
    EXPECT(solver.ordering(), zeroDroppedOrder);
    EXPECT(bitsOf(x), zeroDroppedBits);
}

// Leaves every round to the driver (greedy rounding).
class NoAnswer : public solver::DriverCallback {
public:
    void onRequest(solver::DriverRequest&) override {}
    void onAllFixed() override { ++allFixed; }
    int allFixed = 0;
};

void driverFixing()
{
    QuietProgress q;
    solver::CompactMatrix m = compactOf(gridCase);
    std::vector<double> rhs = gridCase.rhs;
    NoAnswer callback;
    solver::ConstrainedDriver driver(m, rhs, callback, q.progress);
    driver.setIntegerUnknowns({1, 3});
    std::vector<double> x;
    EXPECT(driver.run(x), true);
    EXPECT(driver.status() == solver::ConstrainedDriver::Status::Ok, true);
    EXPECT(callback.allFixed, 1);
    EXPECT(bitsOf(x), greedyFixBits);
    // Column 3 is the unit column, its couplings are zero in both directions.
    bool unit = true;
    for (int k = m.columnStart(3); k < m.columnEnd(3); ++k) {
        unit = unit && m.value(k) == (m.row(k) == 3 ? 1.0 : 0.0);
    }
    EXPECT(unit, true);
}

void driverFixingLists()
{
    QuietProgress q;
    solver::ListMatrix m = listOf(gridCase);
    std::vector<double> rhs = gridCase.rhs;
    NoAnswer callback;
    solver::ConstrainedDriver driver(m, rhs, callback, q.progress);
    driver.setIntegerUnknowns({1, 3});
    std::vector<double> x;
    EXPECT(driver.run(x), true);
    EXPECT(driver.status() == solver::ConstrainedDriver::Status::Ok, true);
    EXPECT(callback.allFixed, 1);
    // Removing the couplings gives the same reduced system as zeroing and dropping them.
    EXPECT(bitsOf(x), greedyFixBits);
    // List 3 is the unit list, no other list refers to unknown 3.
    bool unit = m.list(3).size() == 1 && m.list(3)[0].index == 3 && m.list(3)[0].value == 1.0;
    for (int c = 0; c < m.listCount(); ++c) {
        unit = unit && (c == 3 || solver::ListMatrix::find(m.list(c), 3) < 0);
    }
    EXPECT(unit, true);
}

// Round 0: fixes unknown 0 at its value, stiffens unknowns 2 and 4 and reports
// them as touched; round 1: finished.
class Rounding : public solver::DriverCallback {
public:
    Rounding(solver::ListMatrix& m, std::vector<double>& rhs) : m_(m), rhs_(rhs) {}
    void setDriver(solver::ConstrainedDriver& driver) { driver_ = &driver; }

    void onRequest(solver::DriverRequest& request) override
    {
        request.answer = solver::DriverRequest::Answer::Delivered;
        if (request.round > 0) {
            request.finished = true;
            return;
        }
        driver_->fixVariable(0, (*request.x)[0]);
        for (const auto& [i, target] : {std::pair{2, 1.5}, std::pair{4, -0.5}}) {
            for (solver::SparseEntry& e : m_.list(i)) {
                if (e.index == i) {
                    e.value = e.value + 2.0;
                }
            }
            rhs_[static_cast<std::size_t>(i)] = target * 2.0 + rhs_[static_cast<std::size_t>(i)];
        }
        request.extra = &touched_;
        request.systemChanged = true;
    }

private:
    solver::ListMatrix& m_;
    std::vector<double>& rhs_;
    solver::ConstrainedDriver* driver_ = nullptr;
    std::vector<int> touched_{2, 4, 2};
};

void driverLocalUpdate()
{
    QuietProgress q;
    solver::ListMatrix m = listOf(arrowCase);
    std::vector<double> rhs = arrowCase.rhs;
    Rounding callback(m, rhs);
    solver::DriverOptions options;
    options.keepGoing = true;
    options.localUpdate = true;
    options.local.tolerance = static_cast<double>(0.1f);
    options.local.maxIterations = 50;
    options.fullSolveEvery = 250;
    solver::ConstrainedDriver driver(m, rhs, callback, q.progress, options);
    callback.setDriver(driver);
    std::vector<double> x;
    EXPECT(driver.run(x), true);
    EXPECT(bitsOf(x), cgUpdateBits);
    EXPECT(m.list(0).size(), std::size_t{1});
}

} // namespace

// |A v - value v| of a symmetric matrix [xx, xy, xz, yy, yz, zz].
double eigenResidual(const std::array<double, 6>& m, const field::SymmetricEigen3::Pair& pair)
{
    const Vec3d& v = pair.vector;
    const Vec3d av{m[0] * v.x + m[1] * v.y + m[2] * v.z, m[1] * v.x + m[3] * v.y + m[4] * v.z,
                   m[2] * v.x + m[4] * v.y + m[5] * v.z};
    return length(Vec3d{av.x - pair.value * v.x, av.y - pair.value * v.y, av.z - pair.value * v.z});
}

void jacobiEigen()
{
    using field::SymmetricEigen3;
    // Zero matrix: no rotation, the last pair is (0, 0, 1) with value 0.
    const auto zero = SymmetricEigen3::solve({0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    EXPECT(zero[2].vector.x == 0.0 && zero[2].vector.y == 0.0 && zero[2].vector.z == 1.0 && zero[2].value == 0.0, true);
    // Diagonal: sorted by magnitude, each value keeps its own axis.
    const auto diagonal = SymmetricEigen3::solve({2.0, 0.0, 0.0, -3.0, 0.0, 1.0});
    EXPECT(diagonal[0].value == 1.0 && diagonal[1].value == 2.0 && diagonal[2].value == -3.0, true);
    EXPECT(diagonal[0].vector.z == 1.0 && diagonal[1].vector.x == 1.0 && diagonal[2].vector.y == 1.0, true);
    // Equal magnitudes keep the order xx, yy, zz.
    const auto tie = SymmetricEigen3::solve({-2.0, 0.0, 0.0, 2.0, 0.0, 1.0});
    EXPECT(tie[1].value == -2.0 && tie[1].vector.x == 1.0 && tie[2].value == 2.0 && tie[2].vector.y == 1.0, true);
    // One rotation in the xy plane: the values 1 and 3 come out exactly, 5 stays.
    const std::array<double, 6> plane{2.0, 1.0, 0.0, 2.0, 0.0, 5.0};
    const auto rotated = SymmetricEigen3::solve(plane);
    EXPECT(rotated[0].value == 1.0 && rotated[1].value == 3.0 && rotated[2].value == 5.0, true);
    EXPECT(rotated[2].vector.z == 1.0, true);
    EXPECT(eigenResidual(plane, rotated[0]) < 1e-15 && eigenResidual(plane, rotated[1]) < 1e-15, true);
    // A full matrix over several sweeps: every pair solves A v = value v, the
    // vectors are orthonormal, the magnitudes ascend and the trace is kept.
    const std::array<double, 6> full{4.0, 1.0, -2.0, 2.0, 0.5, -3.0};
    const auto pairs = SymmetricEigen3::solve(full);
    bool solved = true;
    for (int i = 0; i < 3; ++i) {
        solved = solved && eigenResidual(full, pairs[static_cast<std::size_t>(i)]) < 1e-12;
        for (int j = 0; j < 3; ++j) {
            const double d = dot(pairs[static_cast<std::size_t>(i)].vector, pairs[static_cast<std::size_t>(j)].vector);
            solved = solved && std::fabs(d - (i == j ? 1.0 : 0.0)) < 1e-12;
        }
    }
    EXPECT(solved, true);
    EXPECT(std::fabs(pairs[0].value) <= std::fabs(pairs[1].value) &&
               std::fabs(pairs[1].value) <= std::fabs(pairs[2].value), true);
    EXPECT(std::fabs((pairs[0].value + pairs[1].value + pairs[2].value) - 3.0) < 1e-12, true);
}

void curvatureCylinder()
{
    // Open cylinder of radius 1 along z: 24 x 16 squares of two triangles. Only
    // the edges along the axis bend, so the direction of least curvature runs
    // along the axis; a band around its middle fixes the triangles beside it.
    constexpr int ring = 24;
    constexpr int rows = 16;
    std::vector<Vec3f> p;
    for (int j = 0; j <= rows; ++j) {
        for (int i = 0; i < ring; ++i) {
            const double a = 2.0 * std::numbers::pi * i / ring;
            p.push_back({static_cast<float>(std::cos(a)), static_cast<float>(std::sin(a)), 0.25f * static_cast<float>(j)});
        }
    }
    std::vector<int> corners;
    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < ring; ++i) {
            const int a = j * ring + i;
            const int b = j * ring + (i + 1) % ring;
            corners.insert(corners.end(), {a, b, b + ring, a, b + ring, a + ring});
        }
    }
    parts::Part part;
    part.mesh.swapPoints(p);
    part.mesh.editTopology().swapCorners(corners, FaceMode::Triangles);
    part.mesh.topologyChanged();
    part.quadBudget = 400;
    curves::GuideCurve band;
    band.kind = ARM_CURVE_BAND;
    for (int i = 0; i <= 48; ++i) {
        const double a = 2.0 * std::numbers::pi * (i % 48) / 48.0;
        band.points.push_back({static_cast<float>(std::cos(a)), static_cast<float>(std::sin(a)), 2.0f});
    }
    band.source.points = band.points;
    part.curves.push_back(band);
    field::PartEdgeLength::assign(part, 400, 1);
    EXPECT(part.edgeLength > 0.0f, true);
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    field::CurvatureField::compute(part, field::CurvatureField::Options{}, progress);
    EXPECT(part.curvature.has_value(), true);
    if (!part.curvature) {
        return;
    }
    const field::Curvature& c = *part.curvature;
    EXPECT(part.curves[0].width > 0.0f && part.curves[0].cutPoints.size() == 49, true);
    int along = 0;
    int free = 0;
    int fixedAround = 0;
    int fixed = 0;
    for (int f = 0; f < part.mesh.faceCount(); ++f) {
        const std::size_t t = static_cast<std::size_t>(f);
        const float z = part.mesh.points()[static_cast<std::size_t>(part.mesh.topology().vertex(3 * f))].z;
        if ((c.flags[t] & field::curvature_flags::fixedMask) != 0) {
            ++fixed;
            const bool curve = (c.flags[t] & field::curvature_flags::fixedMask) == 4;
            fixedAround += curve && c.weight[t] == 50000.0f && std::fabs(c.direction[t].z) < 0.2f ? 1 : 0;
        } else if (z >= 0.75f && z <= 3.0f && std::fabs(z - 2.0f) > 0.6f) {
            ++free;
            along += (c.flags[t] & field::curvature_flags::direction) != 0 && std::fabs(c.direction[t].z) > 0.99f &&
                             c.weight[t] >= 500.0f ? 1 : 0;
        }
    }
    EXPECT(free > 0 && along == free, true);
    EXPECT(fixed > 0 && fixedAround == fixed, true);
}

// Curvature, cross field and cut graph of a closed triangle mesh. Around every
// closed surface the indices of the singular points add up to its Euler
// characteristic: counted in quarter turns, the sum is 4 (points - edges + faces).
// Returns the part with everything that was computed on it.
parts::Part crossFieldIndexSum(std::vector<Vec3f> points, std::vector<int> corners, int expectedQuarters)
{
    parts::Part part;
    part.mesh.swapPoints(points);
    part.mesh.editTopology().swapCorners(corners, FaceMode::Triangles);
    part.mesh.topologyChanged();
    part.quadBudget = 400;
    field::PartEdgeLength::assign(part, 400, 1);
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    field::CurvatureField::compute(part, field::CurvatureField::Options{}, progress);
    field::CrossFieldSolver::compute(part, field::CrossFieldSolver::Options{}, progress);
    EXPECT(part.crossField.has_value(), true);
    if (!part.crossField) {
        return part;
    }
    const field::CrossField& field = *part.crossField;
    const int euler = part.mesh.pointCount() - part.mesh.adjacency().edgeCount() + part.mesh.faceCount();
    EXPECT(4 * euler, expectedQuarters);
    bool finite = true;
    for (const float angle : field.angle) {
        finite = finite && std::isfinite(angle);
    }
    EXPECT(finite, true);
    int quarters = 0;
    const field::Singularities& singular = field.singularities;
    for (std::size_t i = 0; i < singular.points.size(); ++i) {
        quarters += singular.quarterTurns(i);
    }
    EXPECT(quarters, expectedQuarters);

    // The cut graph makes the part a disc: the cut copy is one piece with
    // points - edges + faces = 1. Across every inner edge that is not cut the
    // combed crosses agree: the quarter turns left there are a multiple of 4.
    field::CutGraphBuilder::compute(part, field::CutGraphBuilder::Options{}, progress);
    EXPECT(part.cutGraph.has_value(), true);
    if (!part.cutGraph) {
        return part;
    }
    const field::CutGraph& cuts = *part.cutGraph;
    const Adjacency& copy = cuts.copy.adjacency();
    EXPECT(cuts.copy.pointCount() - copy.edgeCount() + part.mesh.faceCount(), 1);
    EXPECT(ConnectedParts(cuts.copy.topology(), copy).partCount(), 1);
    const Adjacency& adjacency = part.mesh.adjacency();
    int disagreeing = 0;
    for (int e = 0; e < adjacency.edgeCount(); ++e) {
        const bool cut = cuts.cutEdgeBits.testBounded(static_cast<std::size_t>(e));
        if (!cut && adjacency.edgeCorners(e).size() == 2 && cuts.residualJump[static_cast<std::size_t>(e)] % 4 != 0) {
            ++disagreeing;
        }
    }
    EXPECT(disagreeing, 0);
    return part;
}

// Lock of the border edges on the plane of a kept half. A half that is handed
// in as such without naming an axis locks along all three axes, a named axis
// along that axis only. The half mode that a cutting guide curve causes while
// the symmetry is to be found names no axis and locks nothing.
void borderLockAxes()
{
    using settings::SettingsConversion;
    // Unit square in the plane z = 0: all four border edges lie on the plane
    // of the largest z, one of them on the plane of the largest x.
    parts::Part part;
    std::vector<Vec3f> p{{0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}, {1.0f, 1.0f, 0.0f}, {0.0f, 1.0f, 0.0f}};
    std::vector<int> corners{0, 1, 2, 0, 2, 3};
    part.mesh.swapPoints(p);
    part.mesh.editTopology().swapCorners(corners, FaceMode::Triangles);
    part.mesh.topologyChanged();
    auto locked = [&part](const settings::SymmetrySwitches& symmetry) {
        field::BorderLock::Options options = field::BorderLock::optionsFor(symmetry);
        options.planeMax = {1.0f, 1.0f, 0.0f};
        field::BorderLock::assign(part, options);
        return part.lockedBorderEdges.count();
    };
    arm_settings block;
    settings::Settings::fillDefaults(block);
    block.detect_symmetry = 0;
    block.mirror_x = 0;
    block.mirror_y = 0;
    block.mirror_z = 0;
    block.freeze_half = 1;
    EXPECT(locked(SettingsConversion::symmetry(block)), std::size_t{4});
    block.mirror_x = 1;
    EXPECT(locked(SettingsConversion::symmetry(block)), std::size_t{1});
    block.mirror_x = 0;
    block.freeze_half = 0;
    block.detect_symmetry = 1;
    settings::SymmetrySwitches found = SettingsConversion::symmetry(block);
    EXPECT(locked(found), std::size_t{0});
    // The cut of the symmetry switches the half mode on for a cutting guide curve.
    found.halfMode = true;
    EXPECT(locked(found), std::size_t{0});
    EXPECT(part.lockedBorderEdges.empty(), true);
}

// Size field of a part with cut graph. Without adaption every size is the
// target edge length. A strict quad count gives the sizes room to both sides
// of the target edge length, and the counting adapts them until the number of
// quads they amount to lies between 0.97 and 1.05 of the budget, or gives up
// after four adaptions; the fit that follows moves that number again. The
// band is therefore checked for the part handed in here, whose finished field
// stays in it.
void sizeFieldCount(parts::Part part)
{
    if (!part.cutGraph) {
        return;
    }
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    field::SizeFieldBuilder::Options even;
    even.densityRatio = 1.0f;
    field::SizeFieldBuilder::compute(part, even, progress);
    EXPECT(part.sizeField.has_value(), true);
    if (!part.sizeField) {
        return;
    }
    int different = 0;
    for (std::size_t t = 0; t < part.sizeField->sizeU.size(); ++t) {
        if (part.sizeField->sizeU[t] != part.edgeLength || part.sizeField->sizeV[t] != part.edgeLength) {
            ++different;
        }
    }
    EXPECT(part.sizeField->sizeU.size(), static_cast<std::size_t>(part.mesh.faceCount()));
    EXPECT(different, 0);

    field::SizeFieldBuilder::Options strict;
    strict.exactCount = true;
    field::SizeFieldBuilder::compute(part, strict, progress);
    EXPECT(part.sizeField.has_value(), true);
    if (!part.sizeField) {
        return;
    }
    const float room = std::sqrt(strict.densityRatio);
    EXPECT(part.sizeField->sizeMin, part.edgeLength / room);
    EXPECT(part.sizeField->sizeMax, part.edgeLength * room);
    const float share = static_cast<float>(part.sizeField->estimatedQuads) / static_cast<float>(part.quadBudget);
    EXPECT(share >= 0.97f && share <= 1.05f, true);
}

// Chains of a closed part with cut graph: every cut edge lies in exactly one
// chain, chains end at nodes and no singular point lies inside one, the side
// lists hold copies of the chain points (two different ones inside a chain),
// the quarter turns across a chain are the same along it, and around every
// node the turns of the chains that end there add up to the index of a
// singular point (to nothing at any other node). A part without singular
// point (`torus`) takes its first node in their place and has chains with
// offset variables along its disc cuts.
void cutChains(parts::Part& part, bool torus)
{
    if (!part.cutGraph) {
        return;
    }
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    param::ChainBuilder::compute(part, param::ChainBuilder::Options{}, progress);
    EXPECT(part.chains.has_value(), true);
    if (!part.chains) {
        return;
    }
    const param::Chains& chains = *part.chains;
    const field::CutGraph& cuts = *part.cutGraph;
    const Adjacency& a = part.mesh.adjacency();
    const std::vector<int>& origin = cuts.copy.pointOrigins();

    std::vector<int> timesChained(static_cast<std::size_t>(a.edgeCount()), 0);
    for (const int e : chains.edges.entries) {
        ++timesChained[static_cast<std::size_t>(e)];
    }
    int wronglyChained = 0;
    for (int e = 0; e < a.edgeCount(); ++e) {
        const int expected = cuts.cutEdgeBits.testBounded(static_cast<std::size_t>(e)) ? 1 : 0;
        wronglyChained += timesChained[static_cast<std::size_t>(e)] != expected ? 1 : 0;
    }
    EXPECT(wronglyChained, 0);

    std::vector<int> indexOf(static_cast<std::size_t>(part.mesh.pointCount()), -1);
    for (std::size_t i = 0; i < chains.singularities.points.size(); ++i) {
        indexOf[static_cast<std::size_t>(chains.singularities.points[i])] = chains.singularities.index[i];
    }
    BitArray isNode(static_cast<std::size_t>(part.mesh.pointCount()));
    for (const int n : chains.nodes) {
        isNode.set(static_cast<std::size_t>(n));
    }
    int looseEnds = 0;
    int singularInside = 0;
    int sideMisfits = 0;
    int unevenTurns = 0;
    std::vector<int> turnsAround(static_cast<std::size_t>(part.mesh.pointCount()), 0);
    for (int c = 0; c < chains.count(); ++c) {
        const auto pts = chains.points.row(c);
        const auto edges = chains.edges.row(c);
        const std::size_t ci = static_cast<std::size_t>(c);
        const bool atNodes = isNode.test(static_cast<std::size_t>(pts.front())) &&
                             isNode.test(static_cast<std::size_t>(pts.back()));
        looseEnds += !atNodes && pts.front() != pts.back() ? 1 : 0;
        for (std::size_t k = 1; k + 1 < pts.size(); ++k) {
            singularInside += indexOf[static_cast<std::size_t>(pts[k])] != -1 ? 1 : 0;
        }
        const std::vector<int>& sideA = chains.sideA[ci];
        const std::vector<int>& sideB = chains.sideB[ci];
        const auto pairs = param::ChainSides::corners(a, edges, pts);
        if (!pairs || sideA.size() != pts.size() || sideB.size() != pts.size()) {
            ++sideMisfits;
            continue;
        }
        for (std::size_t k = 0; k < pts.size(); ++k) {
            const bool copies = origin[static_cast<std::size_t>(sideA[k])] == pts[k] &&
                                origin[static_cast<std::size_t>(sideB[k])] == pts[k];
            const bool apart = k == 0 || k + 1 == pts.size() || sideA[k] != sideB[k];
            sideMisfits += copies && apart ? 0 : 1;
        }
        const int turns = static_cast<int>(chains.rotation[ci]);
        // Quarter turns across chain edge k from side A, and clockwise about `point`.
        auto across = [&](std::size_t k) {
            const int jump = cuts.residualJump[static_cast<std::size_t>(edges[k])];
            return pairs->a[k] == a.edgeCorners(edges[k])[0] ? jump : -jump;
        };
        auto about = [&](std::size_t k, int point) { return a.vertex(pairs->a[k]) == point ? across(k) : -across(k); };
        for (std::size_t k = 0; k < edges.size(); ++k) {
            unevenTurns += (across(k) - turns) % 4 != 0 ? 1 : 0;
        }
        turnsAround[static_cast<std::size_t>(pts.front())] += about(0, pts.front());
        turnsAround[static_cast<std::size_t>(pts.back())] += about(edges.size() - 1, pts.back());
    }
    EXPECT(looseEnds, 0);
    EXPECT(singularInside, 0);
    EXPECT(sideMisfits, 0);
    EXPECT(unevenTurns, 0);
    int wrongTurns = 0;
    for (const int n : chains.nodes) {
        const int around = ((turnsAround[static_cast<std::size_t>(n)] % 4) + 4) % 4;
        const int index = std::max(indexOf[static_cast<std::size_t>(n)], 0);
        wrongTurns += around == index || around == (4 - index) % 4 ? 0 : 1;
    }
    EXPECT(wrongTurns, 0);

    EXPECT(chains.substituted, torus);
    EXPECT(chains.genusChains.empty(), !torus);
    if (torus) {
        EXPECT(chains.singularities.points, (std::vector<int>{chains.nodes.front()}));
        EXPECT(chains.singularNodes.count(), std::size_t{1});
        int wrongOffsets = 0;
        for (std::size_t j = 0; j < chains.genusChains.size(); ++j) {
            const int expected = 2 * cuts.copy.pointCount() + 2 * static_cast<int>(j);
            wrongOffsets += chains.offsetVariable[static_cast<std::size_t>(chains.genusChains[j])] != expected ? 1 : 0;
        }
        EXPECT(wrongOffsets, 0);
    } else {
        // Every singular point lies on the cuts and is a node.
        EXPECT(chains.singularNodes.count(), chains.singularities.points.size());
    }
}

// Flat grid of 6 x 5 points (point = 6 * row + column), every square split
// along its rising diagonal.
void flatGrid(Mesh& mesh)
{
    std::vector<Vec3f> p;
    for (int j = 0; j < 5; ++j) {
        for (int i = 0; i < 6; ++i) {
            p.push_back({static_cast<float>(i), static_cast<float>(j), 0.0f});
        }
    }
    std::vector<int> corners;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 5; ++i) {
            const int a = 6 * j + i;
            corners.insert(corners.end(), {a, a + 1, a + 7, a, a + 7, a + 6});
        }
    }
    mesh.swapPoints(p);
    mesh.editTopology().swapCorners(corners, FaceMode::Triangles);
    mesh.topologyChanged();
}

std::vector<int> edgesBetween(const Adjacency& a, std::initializer_list<std::pair<int, int>> pairs)
{
    std::vector<int> edges;
    for (const auto& [p, q] : pairs) {
        edges.push_back(a.edgeBetween(p, q));
    }
    return edges;
}

// First point of every chain.
std::vector<int> chainStarts(const Csr& chains)
{
    std::vector<int> starts;
    for (int r = 0; r < chains.rows(); ++r) {
        starts.push_back(chains.row(r)[0]);
    }
    return starts;
}

// Start rules of the chains that callers can choose: the order of an edge
// list decides which end starts first; starts that follow the chains begin at
// the first branch point instead of the last; the start of pure loops is a
// node once, also when it separates.
void chainStartRules()
{
    Mesh m;
    flatGrid(m);
    const Adjacency& a = m.adjacency();
    const std::vector<int> path = edgesBetween(a, {{7, 13}, {13, 19}});
    const std::vector<int> back(path.rbegin(), path.rend());
    EXPECT(EdgeChains::build(a, path, EdgeChains::Options{}).points.entries, (std::vector<int>{19, 13, 7}));
    EXPECT(EdgeChains::build(a, back, EdgeChains::Options{}).points.entries, (std::vector<int>{7, 13, 19}));
    // An edge listed again counts once, at its first place.
    const std::vector<int> twice{path[1], path[0], path[1]};
    EXPECT(EdgeChains::build(a, twice, EdgeChains::Options{}).points.entries, (std::vector<int>{7, 13, 19}));

    // Three ways between the branch points 7 and 9.
    const std::vector<int> ways = edgesBetween(
        a, {{7, 8}, {8, 9}, {1, 7}, {1, 2}, {2, 3}, {3, 9}, {7, 13}, {13, 14}, {14, 15}, {9, 15}});
    BitArray chosen(static_cast<std::size_t>(a.edgeCount()));
    for (const int e : ways) {
        chosen.set(static_cast<std::size_t>(e));
    }
    EdgeChains::Options following;
    following.startsFollowChains = true;
    EXPECT(chainStarts(EdgeChains::build(a, chosen)), (std::vector<int>{9, 9, 7}));
    const EdgeChains::Result followed = EdgeChains::build(a, chosen, following);
    EXPECT(chainStarts(followed.points), (std::vector<int>{7, 7, 9}));
    EXPECT(followed.nodes, (std::vector<int>{7, 9}));

    // A loop around one square, its first point listed as separator.
    BitArray loop(static_cast<std::size_t>(a.edgeCount()));
    for (const int e : edgesBetween(a, {{7, 8}, {8, 14}, {13, 14}, {7, 13}})) {
        loop.set(static_cast<std::size_t>(e));
    }
    const std::vector<int> separators{7};
    EdgeChains::Options separated;
    separated.separators = separators;
    const EdgeChains::Result around = EdgeChains::build(a, loop, separated);
    EXPECT(around.nodes, (std::vector<int>{7}));
    EXPECT(around.points.starts, (std::vector<int>{0, 5}));
    EXPECT(around.points.entries.front() == 7 && around.points.entries.back() == 7, true);
}

// Chains of an open part with a given cut 21-20-26-19-13: point 26 lies on
// the border, point 20 is singular; both separate. A curve 7-8-9 away from
// the cut runs along the first combed direction, which is turned around in
// the triangles of its second edge, with a curve check cut on either side of
// the curve at point 8: the curve keeps both edges and is split at point 8.
void chainsOnOpenMesh()
{
    parts::Part part;
    flatGrid(part.mesh);
    const Adjacency& a = part.mesh.adjacency();
    const std::size_t faces = static_cast<std::size_t>(part.mesh.faceCount());
    const std::size_t edgeCount = static_cast<std::size_t>(a.edgeCount());

    field::CrossField cross;
    cross.frameE1.assign(faces, Vec3f{1.0f, 0.0f, 0.0f});
    cross.frameE2.assign(faces, Vec3f{0.0f, 1.0f, 0.0f});
    cross.singularities.points = {20};
    cross.singularities.index = {1};
    cross.singularities.fullTurn = {0};

    field::CutGraph cuts;
    cuts.cutEdgeBits = BitArray(edgeCount);
    for (const int e : edgesBetween(a, {{20, 21}, {20, 26}, {19, 26}, {13, 19}})) {
        cuts.cutEdgeBits.set(static_cast<std::size_t>(e));
    }
    cuts.cutEdgeBits.forEachSet([&](std::size_t e) { cuts.cutEdges.push_back(static_cast<int>(e)); });
    cuts.copy = field::CutCopy::cut(part.mesh, cuts.cutEdges, cuts.cutEdgeBits);
    cuts.residualJump.assign(edgeCount, 0);
    cuts.combedU.assign(faces, std::array<float, 2>{1.0f, 0.0f});
    cuts.combedV.assign(faces, std::array<float, 2>{0.0f, 1.0f});
    cuts.curveEdges = edgesBetween(a, {{7, 8}, {8, 9}});
    cuts.curveEdgeKind.assign(edgeCount, field::curve_edge::none);
    for (const int e : cuts.curveEdges) {
        cuts.curveEdgeKind[static_cast<std::size_t>(e)] = field::curve_edge::edgeLoop;
    }
    for (const int c : a.edgeCorners(cuts.curveEdges[1])) {
        cuts.combedU[static_cast<std::size_t>(c / 3)] = {-1.0f, 0.0f};
    }
    cuts.curveCheckCutBits = BitArray(edgeCount);
    for (const int e : edgesBetween(a, {{2, 8}, {8, 14}})) {
        cuts.curveCheckCutBits.set(static_cast<std::size_t>(e));
    }
    field::CurveAlignment::orient(part.mesh, cross, cuts);
    part.crossField = std::move(cross);
    part.cutGraph = std::move(cuts);

    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    param::ChainBuilder::compute(part, param::ChainBuilder::Options{}, progress);
    EXPECT(part.chains.has_value(), true);
    if (!part.chains) {
        return;
    }
    const param::Chains& chains = *part.chains;
    EXPECT(chains.nodes, (std::vector<int>{13, 21, 20, 26}));
    EXPECT(chains.singularNodes.test(2) && chains.singularNodes.count() == 1, true);
    EXPECT(chains.points.starts, (std::vector<int>{0, 2, 5, 7}));
    EXPECT(chains.points.entries, (std::vector<int>{21, 20, 13, 19, 26, 20, 26}));
    const std::vector<int>& origin = part.cutGraph->copy.pointOrigins();
    int sideMisfits = 0;
    for (int c = 0; c < chains.count(); ++c) {
        const auto pts = chains.points.row(c);
        const std::vector<int>& sideA = chains.sideA[static_cast<std::size_t>(c)];
        const std::vector<int>& sideB = chains.sideB[static_cast<std::size_t>(c)];
        if (sideA.size() != pts.size() || sideB.size() != pts.size()) {
            ++sideMisfits;
            continue;
        }
        for (std::size_t k = 0; k < pts.size(); ++k) {
            sideMisfits += origin[static_cast<std::size_t>(sideA[k])] == pts[k] &&
                                   origin[static_cast<std::size_t>(sideB[k])] == pts[k]
                               ? 0
                               : 1;
        }
    }
    EXPECT(sideMisfits, 0);
    // Inside the chain 13-19-26 the cut separates the two copies of point 19.
    EXPECT(chains.sideA[1][1] != chains.sideB[1][1], true);

    const param::CurveChains& curves = chains.curves;
    EXPECT(std::count(curves.edgeKind.begin(), curves.edgeKind.end(), field::curve_edge::edgeLoop), 2);
    EXPECT(curves.edges.starts, (std::vector<int>{0, 1, 2}));
    EXPECT(curves.coordinate, (std::vector<std::int8_t>{1, 1}));
    EXPECT(curves.turnPoints, (std::vector<int>{8}));
}

// Curve chains on the flat grid without cuts. The frames and the combed
// directions lie along the axes; in the triangles `turned` the first combed
// direction is `first` and the second one stands a quarter turn left of it.
// `checkCuts` are curve check cuts.
param::CurveChains curveChainsOnGrid(const Mesh& mesh, field::CutGraph& cuts, const std::vector<int>& curveEdges,
                                     const std::vector<int>& turned, std::array<float, 2> first,
                                     const std::vector<int>& checkCuts)
{
    const Adjacency& a = mesh.adjacency();
    const std::size_t faces = static_cast<std::size_t>(mesh.faceCount());
    const std::size_t edgeCount = static_cast<std::size_t>(a.edgeCount());
    field::CrossField cross;
    cross.frameE1.assign(faces, Vec3f{1.0f, 0.0f, 0.0f});
    cross.frameE2.assign(faces, Vec3f{0.0f, 1.0f, 0.0f});
    cuts.cutEdgeBits = BitArray(edgeCount);
    cuts.copy = field::CutCopy::cut(mesh, cuts.cutEdges, cuts.cutEdgeBits);
    cuts.residualJump.assign(edgeCount, 0);
    cuts.combedU.assign(faces, std::array<float, 2>{1.0f, 0.0f});
    cuts.combedV.assign(faces, std::array<float, 2>{0.0f, 1.0f});
    for (const int t : turned) {
        cuts.combedU[static_cast<std::size_t>(t)] = first;
        cuts.combedV[static_cast<std::size_t>(t)] = {-first[1], first[0]};
    }
    cuts.curveEdges = curveEdges;
    cuts.curveEdgeKind.assign(edgeCount, field::curve_edge::none);
    for (const int e : curveEdges) {
        cuts.curveEdgeKind[static_cast<std::size_t>(e)] = field::curve_edge::edgeLoop;
    }
    cuts.curveCheckCutBits = BitArray(edgeCount);
    for (const int e : checkCuts) {
        cuts.curveCheckCutBits.set(static_cast<std::size_t>(e));
    }
    field::CurveAlignment::orient(mesh, cross, cuts);
    return param::CurveChainBuilder::build(mesh, cross, cuts, false);
}

// The one edge of the cut copy over edge `edge` of a mesh without cuts.
int copyEdge(const Mesh& mesh, const field::CutGraph& cuts, int edge)
{
    const std::vector<int> copies = cuts.copy.copyEdges(mesh.adjacency(), edge);
    EXPECT(copies.size(), std::size_t{1});
    return copies.empty() ? -1 : copies.front();
}

// Two curve edges leave point 8 to the right (to 9) and upwards (to 14). In
// the two triangles of the edge 8-14 the first combed direction points up,
// so both edges run along the first direction, away from point 8: label 0.
// Around point 8 the walk reaches one of them with the label of the other
// across the quarter sector between them (one inner edge, 8-15, which does
// not point away); the other way round the sector is three quarters wide
// and bends back. Whichever of the two the walk starts at, it comes to the
// edge 8-14 over the quarter sector: that edge gets no chain. At the points
// 9 and 14 the walk only comes back to its own edge around the whole point.
void curveChainExclusion()
{
    Mesh mesh;
    flatGrid(mesh);
    const Adjacency& a = mesh.adjacency();
    const std::vector<int> curve = edgesBetween(a, {{8, 9}, {8, 14}});
    std::vector<int> turned;
    for (const int c : a.edgeCorners(curve[1])) {
        turned.push_back(c / 3);
    }
    field::CutGraph cuts;
    const param::CurveChains chains = curveChainsOnGrid(mesh, cuts, curve, turned, {0.0f, 1.0f}, {});
    const int kept = copyEdge(mesh, cuts, curve[0]);
    const int excluded = copyEdge(mesh, cuts, curve[1]);
    EXPECT(chains.edges.starts, (std::vector<int>{0, 1}));
    EXPECT(chains.edges.entries, (std::vector<int>{kept}));
    EXPECT(chains.coordinate, (std::vector<std::int8_t>{1}));
    EXPECT(chains.turnPoints.empty(), true);
    EXPECT(chains.edgeKind[static_cast<std::size_t>(kept)], field::curve_edge::edgeLoop);
    EXPECT(chains.edgeKind[static_cast<std::size_t>(excluded)], field::curve_edge::none);
}

// A closed run of curve edges around the square 7-8-14-13 and one curve edge
// 10-16 away from it. Left of the edges 8-14 and 7-13 lies the triangle of
// their first corner; there the first combed direction points down, so all
// four edges of the square follow the first direction and the edge 10-16 the
// second. The diagonal 7-14 is a curve check cut, so no edge is excluded at
// its ends; at the points 8 and 13 the labels of the two curve edges differ.
// The edge 10-16 is chained first (second direction). The run starts at
// point 7, the smaller end of its first edge, with the edge 7-8 and goes
// 7 -> 8 -> 14 -> 13 -> 7. Its labels, each read at the first corner of the
// edge: 7-8 along the first direction (0), 8-14 upwards against it (2),
// 14-13 to the left against it (2), 13-7 downwards along it (0). It is split
// at point 8 and at point 13; the last piece (7-13) has the label of the
// first (7-8) and goes in front of it, and the piece behind moves on by one
// edge.
void curveChainClosedRun()
{
    Mesh mesh;
    flatGrid(mesh);
    const Adjacency& a = mesh.adjacency();
    const std::vector<int> curve = edgesBetween(a, {{7, 8}, {8, 14}, {13, 14}, {7, 13}, {10, 16}});
    const std::vector<int> turned{a.edgeCorners(curve[1])[0] / 3, a.edgeCorners(curve[3])[0] / 3};
    field::CutGraph cuts;
    const param::CurveChains chains =
        curveChainsOnGrid(mesh, cuts, curve, turned, {0.0f, -1.0f}, edgesBetween(a, {{7, 14}}));
    std::vector<int> onCopy;
    for (const int e : curve) {
        onCopy.push_back(copyEdge(mesh, cuts, e));
    }
    EXPECT(chains.edges.starts, (std::vector<int>{0, 1, 3, 5}));
    EXPECT(chains.edges.entries, (std::vector<int>{onCopy[4], onCopy[3], onCopy[0], onCopy[1], onCopy[2]}));
    EXPECT(chains.coordinate, (std::vector<std::int8_t>{0, 1, 1}));
    EXPECT(chains.turnPoints, (std::vector<int>{8, 13}));
}

namespace {
void systemOnClosedPart(parts::Part& part, bool torus);
void solutionOnClosedPart(parts::Part& part);
}

void crossFieldSphere()
{
    // Sphere of radius 1: two poles and 11 rings of 24 points.
    constexpr int ring = 24;
    constexpr int rings = 11;
    std::vector<Vec3f> p{{0.0f, 0.0f, 1.0f}};
    for (int j = 1; j <= rings; ++j) {
        const double polar = std::numbers::pi * j / (rings + 1);
        for (int i = 0; i < ring; ++i) {
            const double a = 2.0 * std::numbers::pi * i / ring;
            p.push_back({static_cast<float>(std::sin(polar) * std::cos(a)), static_cast<float>(std::sin(polar) * std::sin(a)),
                         static_cast<float>(std::cos(polar))});
        }
    }
    p.push_back({0.0f, 0.0f, -1.0f});
    const int south = static_cast<int>(p.size()) - 1;
    std::vector<int> corners;
    for (int i = 0; i < ring; ++i) {
        const int next = (i + 1) % ring;
        corners.insert(corners.end(), {0, 1 + i, 1 + next});
        const int last = 1 + (rings - 1) * ring;
        corners.insert(corners.end(), {south, last + next, last + i});
        for (int j = 0; j + 1 < rings; ++j) {
            const int a = 1 + j * ring + i;
            const int b = 1 + j * ring + next;
            corners.insert(corners.end(), {a, a + ring, b + ring, a, b + ring, b});
        }
    }
    parts::Part part = crossFieldIndexSum(std::move(p), std::move(corners), 8);
    cutChains(part, false);
    systemOnClosedPart(part, false);
    solutionOnClosedPart(part);
    sizeFieldCount(std::move(part));
}

void crossFieldTorus()
{
    // Torus around z: 32 x 16 squares of two triangles, radii 2 and 0.7.
    constexpr int around = 32;
    constexpr int tube = 16;
    std::vector<Vec3f> p;
    for (int j = 0; j < around; ++j) {
        const double a = 2.0 * std::numbers::pi * j / around;
        for (int i = 0; i < tube; ++i) {
            const double b = 2.0 * std::numbers::pi * i / tube;
            const double r = 2.0 + 0.7 * std::cos(b);
            p.push_back({static_cast<float>(r * std::cos(a)), static_cast<float>(r * std::sin(a)),
                         static_cast<float>(0.7 * std::sin(b))});
        }
    }
    std::vector<int> corners;
    for (int j = 0; j < around; ++j) {
        for (int i = 0; i < tube; ++i) {
            const int a = j * tube + i;
            const int b = j * tube + (i + 1) % tube;
            const int c = ((j + 1) % around) * tube + i;
            const int d = ((j + 1) % around) * tube + (i + 1) % tube;
            corners.insert(corners.end(), {a, c, d, a, d, b});
        }
    }
    parts::Part part = crossFieldIndexSum(std::move(p), std::move(corners), 0);
    cutChains(part, true);
    systemOnClosedPart(part, true);
    solutionOnClosedPart(part);
}

namespace {

// Coefficients of the checks that are not finite, and one next to which two
// neighbouring numbers lie 16 apart.
constexpr double notNumber = std::numeric_limits<double>::quiet_NaN();
constexpr double infinite = std::numeric_limits<double>::infinity();
constexpr double big = 72057594037927936.0; // 2^56

// Term lists, rows and conditions flattened to numbers and coefficient bits.
// Every value that is not a number counts as the same one.
void flatten(const param::TermList& terms, std::vector<std::uint64_t>& out)
{
    out.push_back(terms.size());
    for (const param::Term& t : terms) {
        out.push_back(static_cast<std::uint64_t>(static_cast<std::int64_t>(t.index)));
        out.push_back(std::bit_cast<std::uint64_t>(std::isnan(t.value) ? notNumber : t.value));
    }
}

std::vector<std::uint64_t> flat(const param::TermList& terms)
{
    std::vector<std::uint64_t> out;
    flatten(terms, out);
    return out;
}

std::vector<std::uint64_t> flat(const std::vector<param::TermList>& rows)
{
    std::vector<std::uint64_t> out;
    for (const param::TermList& row : rows) {
        flatten(row, out);
    }
    return out;
}

std::vector<std::uint64_t> flat(const std::vector<param::Constraint>& conditions)
{
    std::vector<std::uint64_t> out;
    for (const param::Constraint& c : conditions) {
        out.push_back(static_cast<std::uint64_t>(static_cast<std::int64_t>(c.pivot)));
        flatten(c.terms, out);
    }
    return out;
}

// Column order: variables of the columns, the split and the rows rewritten to columns.
void eliminationColumnOrder()
{
    const std::vector<std::uint8_t> flags{0, 1, 2, 0, 1, 2, 0, 0};
    const std::vector<int> preferred{6, 4};
    const std::vector<param::TermList> rows{{{5, 1.0}, {1, 2.0}, {3, 1.0}}, {{2, 1.0}, {0, 1.0}, {1, 1.0}}};
    param::EliminationSetup setup;
    setup.preferred = preferred;
    setup.flags = flags;
    {
        // Preferred 6, 4; then the variables without the kept-free flag as they
        // occur: 1, 3, 0; behind the split those with it: 5, 2.
        std::vector<param::TermList> work = rows;
        setup.mode = param::EliminationMode::ChainEnds;
        const param::EliminationColumns columns(work, 8, setup);
        EXPECT(v(columns.variables()), (std::vector<int>{6, 4, 1, 3, 0, 5, 2}));
        EXPECT(columns.split(), 5);
        EXPECT(flat(work), flat(std::vector<param::TermList>{{{5, 1.0}, {2, 2.0}, {3, 1.0}},
                                                               {{6, 1.0}, {4, 1.0}, {2, 1.0}}}));
    }
    {
        // The continuous variables of the rows come first: 3, 0; then as before.
        std::vector<param::TermList> work = rows;
        setup.mode = param::EliminationMode::IsoLines;
        const param::EliminationColumns columns(work, 8, setup);
        EXPECT(v(columns.variables()), (std::vector<int>{3, 0, 6, 4, 1, 5, 2}));
        EXPECT(columns.split(), 5);
        EXPECT(flat(work), flat(std::vector<param::TermList>{{{5, 1.0}, {4, 2.0}, {0, 1.0}},
                                                               {{6, 1.0}, {1, 1.0}, {4, 1.0}}}));
    }
    {
        // Preferred 3; variables below the constant column 5: 2, 1; above it:
        // 7, 8; behind the split the constant column itself.
        const std::vector<int> softPreferred{3};
        std::vector<param::TermList> work{{{7, 1.0}, {2, 1.0}, {5, 4.0}}, {{1, 2.0}, {8, 1.0}, {2, 3.0}}};
        param::EliminationSetup soft;
        soft.mode = param::EliminationMode::SoftSystem;
        soft.preferred = softPreferred;
        soft.constantColumn = 5;
        const param::EliminationColumns columns(work, 9, soft);
        EXPECT(v(columns.variables()), (std::vector<int>{3, 2, 1, 7, 8, 5}));
        EXPECT(columns.split(), 5);
        EXPECT(flat(work), flat(std::vector<param::TermList>{{{3, 1.0}, {1, 1.0}, {5, 4.0}},
                                                               {{2, 2.0}, {4, 1.0}, {1, 3.0}}}));
    }
    {
        // Constant column 1. The passes rewrite in place: variable 3 gets
        // column 0 and variable 4 column 1, the number of the constant
        // column. The last pass takes that entry for the constant: both end
        // up in column 2, which stands for variable 1.
        std::vector<param::TermList> work{{{3, 1.0}, {4, 2.0}, {1, 5.0}}};
        param::EliminationSetup soft;
        soft.mode = param::EliminationMode::SoftSystem;
        soft.constantColumn = 1;
        const param::EliminationColumns columns(work, 5, soft);
        EXPECT(v(columns.variables()), (std::vector<int>{3, 4, 1}));
        EXPECT(columns.split(), 2);
        EXPECT(flat(work), flat(std::vector<param::TermList>{{{0, 1.0}, {2, 2.0}, {2, 5.0}}}));
    }
}

// One elimination with its outcome. The expected values are computed step by
// step from the rules, each system is small enough to follow by hand.
struct EliminationCase {
    param::EliminationMode mode;
    int variables;
    std::vector<std::uint8_t> flags;
    std::vector<int> preferred;
    int constantColumn;
    std::vector<param::TermList> rows;
    bool ok;
    int firstSmallest;
    std::vector<param::Constraint> conditions;
    std::vector<param::TermList> rowsLeft;
};

// Columns x0, x1, x3, x2. The last pivot (2 in column 2) is divided out
// first; the pivot 5e-11 of column 0 is below the limit. The rows stay,
// the conditions are in column numbers, the unreached one is empty.
const EliminationCase tooSmallPivot{param::EliminationMode::IsoLines, 4, {1, 1, 1, 1}, {0}, -1,
                                    {{{0, 5e-11}, {1, 2.0}}, {{3, 2.0}, {2, 3.0}}},
                                    false, 0, {{-1, {}}, {2, {{3, -1.5}}}},
                                    {{{0, 5e-11}, {1, 2.0}}, {{2, 2.0}, {3, 3.0}}}};

const EliminationCase eliminationCases[] = {
    // Leading 1 in a later row: row 1 becomes the pivot row of column 0 and
    // row 0 turns into 3 x1 + x2 - 2 x3. Its 1 stands in the column of x2,
    // which is exchanged with the column of x1: x2 = -3 x1 + 2 x3, x0 = x1 - x3.
    {param::EliminationMode::ChainEnds, 4, {1, 1, 1, 1}, {}, -1,
     {{{0, 2.0}, {1, 1.0}, {2, 1.0}}, {{0, 1.0}, {1, -1.0}, {3, 1.0}}},
     true, -1, {{0, {{1, 1.0}, {3, -1.0}}}, {2, {{1, -3.0}, {3, 2.0}}}},
     {}},
    // Back substitution: x1 = -3 x2 is put into x0 = -2 x1 - x3. The terms
    // keep the order of the columns (x3 got its column before x2).
    {param::EliminationMode::ChainEnds, 4, {1, 1, 1, 1}, {}, -1,
     {{{0, 1.0}, {1, 2.0}, {3, 1.0}}, {{1, 1.0}, {2, 3.0}}},
     true, -1, {{0, {{3, -1.0}, {2, 6.0}}}, {1, {{2, -3.0}}}},
     {}},
    // No coefficient of magnitude 1: column 0 takes the smallest magnitude, 2
    // in rows 1 and 2, the first of them wins. Row 0 becomes -12 x1 - 6 x3 -
    // 12 x2 and is divided by 6; row 2 becomes -10 x1 - 14 x2 and stays (14 / 10
    // is not whole). The -1 of x3 is exchanged into column 1; column 2 takes
    // the -10 of x1. Solution x1 = -1.4 x2, x3 = (2.8 - 2) x2, x0 = 1.5 x2.
    {param::EliminationMode::ChainEnds, 4, {1, 1, 1, 1}, {}, -1,
     {{{0, 4.0}, {1, 6.0}, {3, 3.0}}, {{0, -2.0}, {2, 3.0}}, {{0, 2.0}, {1, 5.0}, {2, 4.0}}},
     true, 0, {{0, {{2, 1.5}}}, {3, {{2, 0.7999999999999998}}}, {1, {{2, -1.4}}}},
     {}},
    // Columns x0, x2, x1; only x0 lies in the continuous leading part. The 1
    // of x2 is out of reach for the continuous column, so 2 is the pivot.
    {param::EliminationMode::ChainEnds, 3, {0, 0, 1}, {}, -1,
     {{{0, 2.0}, {2, 1.0}, {1, 3.0}}},
     true, 0, {{0, {{2, -0.5}, {1, -1.5}}}},
     {}},
    // The same row with an integer x0: its column may be exchanged up to the
    // split, the 1 of x2 becomes the pivot.
    {param::EliminationMode::ChainEnds, 3, {1, 0, 1}, {}, -1,
     {{{0, 2.0}, {2, 1.0}, {1, 3.0}}},
     true, -1, {{2, {{0, -2.0}, {1, -3.0}}}},
     {}},
    // Every column continuous: the leading part is all columns, the 1 of x1
    // is exchanged into column 0.
    {param::EliminationMode::ChainEnds, 3, {0, 0, 0}, {}, -1,
     {{{0, 2.0}, {1, 1.0}, {2, 3.0}}},
     true, -1, {{1, {{0, -2.0}, {2, -3.0}}}},
     {}},
    // Column 0 belongs to the preferred variable 3, which occurs in no row:
    // no pivot there, the row stays for column 1.
    {param::EliminationMode::ChainEnds, 4, {1, 1, 1, 1}, {3}, -1,
     {{{0, 2.0}, {1, 3.0}}},
     true, 0, {{0, {{1, -1.5}}}},
     {}},
    // x1 has the kept-free flag: its column lies behind the split, its 1 is
    // not exchanged into column 0.
    {param::EliminationMode::ChainEnds, 3, {1, 2, 1}, {}, -1,
     {{{0, 2.0}, {1, 1.0}, {2, 4.0}}},
     true, 0, {{0, {{2, -2.0}, {1, -0.5}}}},
     {}},
    // Columns x1, x3 (continuous), x2 (preferred), x0. Column 0 finds the 1 of
    // x3 in row 1 within the continuous part; column 1 then has only the 2 of
    // x1, the pivot of condition 1.
    {param::EliminationMode::IsoLines, 4, {1, 0, 1, 0}, {2}, -1,
     {{{0, 1.0}, {1, 2.0}, {2, 1.0}}, {{3, 1.0}, {0, 1.0}}},
     true, 1, {{3, {{0, -1.0}}}, {1, {{2, -0.5}, {0, -0.5}}}},
     {}},
    tooSmallPivot,
    // Constant column 6. Row 0 leads with 1. Row 1 becomes -4 x1 - x2 + 4 x7 -
    // 6 c, its -1 is exchanged into column 1. Row 2 holds only the constant
    // and is passed over. Row 3 has no 1: its first coefficient 0.5 is the pivot.
    {param::EliminationMode::SoftSystem, 8, {}, {}, 6,
     {{{0, 1.0}, {1, 2.0}, {6, 3.0}}, {{0, 2.0}, {2, -1.0}, {7, 4.0}}, {{6, 2.0}}, {{1, 0.5}, {3, 4.0}}},
     true, -1, {{0, {{3, 16.0}, {6, -3.0}}}, {2, {{3, 32.0}, {7, 4.0}, {6, -6.0}}}, {1, {{3, -8.0}}}},
     {}},
    // The pivot 5e-11 of row 0 is below the limit: its condition is left out.
    // In x1 = -1e-9 x2 - 2 c the tiny term is dropped.
    {param::EliminationMode::SoftSystem, 4, {}, {}, 3,
     {{{0, 5e-11}, {1, 2.0}}, {{1, 1.0}, {2, 1e-09}, {3, 2.0}}},
     true, -1, {{1, {{3, -2.0}}}},
     {}},
    // Columns x0 (continuous), x1, x2. Column 0 gets no pivot: its only
    // coefficient is 0. For column 1 the 1 of x2 is exchanged into it, but the
    // row still leads with the 0 of column 0, and the leading coefficient is
    // the pivot: below the limit, the rows stay in column numbers.
    {param::EliminationMode::ChainEnds, 3, {0, 1, 1}, {}, -1,
     {{{0, 0.0}, {1, 2.0}, {2, 1.0}}},
     false, -1, {{-1, {}}},
     {{{0, 0.0}, {1, 1.0}, {2, 2.0}}}},
    // A coefficient that is not a number never is the smallest magnitude: the
    // 4 of row 1 is the pivot of column 0. Subtracting leaves row 0 with that
    // coefficient in column 0, so it leads in no later column and gives no
    // condition.
    {param::EliminationMode::ChainEnds, 2, {1, 1}, {}, -1,
     {{{0, notNumber}, {1, 2.0}}, {{0, 4.0}, {1, 6.0}}},
     true, 0, {{0, {{1, -1.5}}}},
     {}},
    // Constant column 2. A coefficient that is not a number is not 0: it is
    // the pivot of the row and not below the limit, and so are the terms.
    {param::EliminationMode::SoftSystem, 3, {}, {}, 2,
     {{{0, notNumber}, {1, 2.0}, {2, 1.0}}},
     true, -1, {{0, {{1, notNumber}, {2, notNumber}}}},
     {}},
};

param::EliminationSetup setupOf(const EliminationCase& c)
{
    param::EliminationSetup setup;
    setup.mode = c.mode;
    setup.preferred = c.preferred;
    setup.flags = c.flags;
    setup.constantColumn = c.constantColumn;
    return setup;
}

void eliminationRules()
{
    QuietProgress q;
    for (const EliminationCase& c : eliminationCases) {
        std::vector<param::TermList> rows = c.rows;
        param::Elimination elimination(rows, c.variables, setupOf(c), q.progress);
        EXPECT(elimination.run(), c.ok);
        EXPECT(elimination.firstSmallestPivot(), c.firstSmallest);
        EXPECT(flat(elimination.constraints()), flat(c.conditions));
        EXPECT(flat(rows), flat(c.rowsLeft));
    }
    // After the failed run only the reached condition is entered into the map
    // variable -> condition.
    const EliminationCase& failed = tooSmallPivot;
    std::vector<param::TermList> rows = failed.rows;
    param::Elimination elimination(rows, failed.variables, setupOf(failed), q.progress);
    EXPECT(elimination.run(), false);
    std::vector<int> conditionOf(4, -1);
    elimination.enterPivots(conditionOf, 10);
    EXPECT(conditionOf, (std::vector<int>{-1, -1, 11, -1}));
    // Asked to, the failed run writes the variables back: columns 2 and 3
    // stand for x3 and x2. The rows and the reached condition x3 = -1.5 x2
    // then speak of variables, their terms in the order of the columns, and
    // the pivot is entered for x3. Asking again changes nothing.
    const std::vector<param::Constraint> inVariables{{-1, {}}, {3, {{2, -1.5}}}};
    const std::vector<param::TermList> rowsInVariables{{{0, 5e-11}, {1, 2.0}}, {{3, 2.0}, {2, 3.0}}};
    for (int pass = 0; pass < 2; ++pass) {
        elimination.restoreVariables();
        EXPECT(flat(elimination.constraints()), flat(inVariables));
        EXPECT(flat(rows), flat(rowsInVariables));
    }
    std::vector<int> restored(4, -1);
    elimination.enterPivots(restored, 10);
    EXPECT(restored, (std::vector<int>{-1, -1, -1, 11}));
}

// x0 - x1 = 0 with `preferredCount` preferred variables in front: the 1 of x0
// is exchanged into column 0 and is the only pivot. With `soft` it is a row
// of the additional system, x0 - x1 + 2 c = 0, with the constants in the
// variable behind the preferred ones; their column does not count.
std::vector<param::Constraint> eliminatedWithPreferred(int preferredCount, bool soft)
{
    QuietProgress q;
    std::vector<int> preferred;
    for (int i = 0; i < preferredCount; ++i) {
        preferred.push_back(2 + i);
    }
    const int constant = 2 + preferredCount;
    const std::vector<std::uint8_t> flags(static_cast<std::size_t>(constant), 1);
    std::vector<param::TermList> rows{{{0, 1.0}, {1, -1.0}}};
    param::EliminationSetup setup;
    setup.preferred = preferred;
    if (soft) {
        rows.front().push_back({constant, 2.0});
        setup.mode = param::EliminationMode::SoftSystem;
        setup.constantColumn = constant;
    } else {
        setup.flags = flags;
    }
    param::Elimination elimination(rows, constant + (soft ? 1 : 0), setup, q.progress);
    EXPECT(elimination.run(), true);
    return elimination.takeConstraints();
}

void eliminationColumnBlocks()
{
    const std::vector<param::Constraint> one{{0, {{1, 1.0}}}};
    const std::vector<param::Constraint> twice{{0, {{1, 1.0}}}, {0, {{1, 1.0}}}};
    // 1000 columns are one block: the run ends with the pivot.
    EXPECT(flat(eliminatedWithPreferred(998, false)), flat(one));
    // 1002 columns: the pivot ends the first block only, the second block
    // takes the row once more.
    EXPECT(flat(eliminatedWithPreferred(1000, false)), flat(twice));
    // The same in the additional system, the column of the constants not counted.
    const std::vector<param::Constraint> softOne{{0, {{1, 1.0}, {1000, -2.0}}}};
    const std::vector<param::Constraint> softTwice{{0, {{1, 1.0}, {1002, -2.0}}}, {0, {{1, 1.0}, {1002, -2.0}}}};
    EXPECT(flat(eliminatedWithPreferred(998, true)), flat(softOne));
    EXPECT(flat(eliminatedWithPreferred(1000, true)), flat(softTwice));
}

// 23 terms of 5 variables (term i belongs to variable 7 i mod 5), the first
// of them 0. More than 16 terms are not sorted stably: the terms of a
// variable are added up in the order 20, 5, 15, 10 (x0), 3, 18, 8, 13 (x1),
// 1, 21, 16, 11, 6 (x2), 14, 9, 19, 4 (x3) and 12, 7, 17, 2, 22 (x4). Next to
// `big` every partial sum is rounded to a multiple of 16, so the sum depends
// on the order: x1 is ((big + 9) - big) + 3 = 19, in the order of the list it
// would be ((big - big) + 3) + 9 = 12. For every variable only the order
// above, or the same with its first two terms exchanged, gives these sums
// (19, 19, 48, 19, 41).
const param::TermList manyTerms{
    {0, 0.0}, {2, 3.0},  {4, -big}, {1, big},  {3, 3.0}, {0, big},  {2, -big}, {4, 3.0},
    {1, -big}, {3, big}, {0, 3.0},  {2, 9.0},  {4, 21.0}, {1, 3.0}, {3, 9.0},  {0, -big},
    {2, big},  {4, big}, {1, 9.0},  {3, -big}, {0, 9.0},  {2, 21.0}, {4, 9.0}};
const std::vector<std::uint64_t> manyTermsBits{0x4033000000000000ull, 0x4033000000000000ull, 0x4048000000000000ull,
                                               0x4033000000000000ull, 0x4044800000000000ull};

void constraintTermRules()
{
    using param::ConstraintTerms;
    using param::TermList;
    // 2 + 0.5 for x1, the zero of x2 is skipped, 1 - 1 for x3 vanishes.
    TermList terms{{3, 1.0}, {1, 2.0}, {3, -1.0}, {2, 0.0}, {1, 0.5}};
    ConstraintTerms::normalize(terms);
    EXPECT(flat(terms), flat(TermList{{1, 2.5}}));
    // A list built from two pairs of equal variables with opposite signs is empty.
    terms = {{14, 1.0}, {14, -1.0}, {10, 1.0}, {10, -1.0}};
    ConstraintTerms::normalize(terms);
    EXPECT(terms.empty(), true);
    terms = manyTerms;
    ConstraintTerms::normalize(terms);
    std::vector<int> variables;
    std::vector<double> sums;
    for (const param::Term& t : terms) {
        variables.push_back(t.index);
        sums.push_back(t.value);
    }
    EXPECT(variables, (std::vector<int>{0, 1, 2, 3, 4}));
    EXPECT(bitsOf(sums), manyTermsBits);
    // A sum that is not a number is not 0 and stays: x2 holds such a term,
    // for x3 infinity meets minus infinity. 1 - 1 for x1 vanishes.
    terms = {{2, notNumber}, {1, 1.0}, {2, 1.0}, {1, -1.0}, {3, infinite}, {3, -infinite}};
    ConstraintTerms::normalize(terms);
    EXPECT(flat(terms), flat(TermList{{2, notNumber}, {3, notNumber}}));

    // 1e-9 and -5e-9 are below 1e-8 and go, together with the zero; 1e-8 stays.
    terms = {{0, 1e-09}, {1, 0.0}, {2, -5e-09}, {3, 1e-08}, {4, 2.0}};
    ConstraintTerms::dropTiny(terms, 1e-8);
    EXPECT(flat(terms), flat(TermList{{3, 1e-08}, {4, 2.0}}));
    // Without a tiny term the zero stays.
    terms = {{1, 0.0}, {4, 2.0}};
    ConstraintTerms::dropTiny(terms, 1e-8);
    EXPECT(flat(terms), flat(TermList{{1, 0.0}, {4, 2.0}}));
    // A coefficient that is not a number is not tiny and no zero.
    terms = {{0, notNumber}, {1, 1e-09}, {2, infinite}, {3, 0.0}};
    ConstraintTerms::dropTiny(terms, 1e-8);
    EXPECT(flat(terms), flat(TermList{{0, notNumber}, {2, infinite}}));

    // -2 (x0 - x1 + x3) + (2 x0 + x1 + x2): x0 vanishes.
    EXPECT(flat(ConstraintTerms::combine(-2.0, {{0, 1.0}, {1, -1.0}, {3, 1.0}}, 1.0, {{0, 2.0}, {1, 1.0}, {2, 1.0}})),
           flat(TermList{{1, 3.0}, {2, 1.0}, {3, -2.0}}));
    // A variable of one list only stays even with a product of 0.
    EXPECT(flat(ConstraintTerms::combine(0.0, {{0, 1.0}}, 1.0, {{1, 2.0}})), flat(TermList{{0, 0.0}, {1, 2.0}}));
    // Sums that are not numbers stay (x0; infinity minus infinity for x1);
    // 10 * 1e308 is beyond the largest number.
    EXPECT(flat(ConstraintTerms::combine(1.0, {{0, notNumber}, {1, infinite}}, 10.0,
                                         {{0, 1.0}, {1, -infinite}, {2, 1e308}})),
           flat(TermList{{0, notNumber}, {1, notNumber}, {2, infinite}}));

    // Divided by 6; not divided (14 / 10); not divided (holds a 1); divided by 0.5.
    terms = {{1, -12.0}, {2, -12.0}, {3, -6.0}};
    ConstraintTerms::reduceToIntegers(terms);
    EXPECT(flat(terms), flat(TermList{{1, -2.0}, {2, -2.0}, {3, -1.0}}));
    terms = {{1, -10.0}, {2, -14.0}};
    ConstraintTerms::reduceToIntegers(terms);
    EXPECT(flat(terms), flat(TermList{{1, -10.0}, {2, -14.0}}));
    terms = {{0, 3.0}, {1, 1.0}};
    ConstraintTerms::reduceToIntegers(terms);
    EXPECT(flat(terms), flat(TermList{{0, 3.0}, {1, 1.0}}));
    terms = {{0, 0.5}, {1, 1.5}};
    ConstraintTerms::reduceToIntegers(terms);
    EXPECT(flat(terms), flat(TermList{{0, 1.0}, {1, 3.0}}));
    // Not divided: a quotient that is not a number, an infinite one, and
    // 3e9, which does not fit an int. 2e9 fits.
    terms = {{0, 0.5}, {1, notNumber}};
    ConstraintTerms::reduceToIntegers(terms);
    EXPECT(flat(terms), flat(TermList{{0, 0.5}, {1, notNumber}}));
    terms = {{0, 0.5}, {1, infinite}};
    ConstraintTerms::reduceToIntegers(terms);
    EXPECT(flat(terms), flat(TermList{{0, 0.5}, {1, infinite}}));
    terms = {{0, 0.5}, {1, 1500000000.0}};
    ConstraintTerms::reduceToIntegers(terms);
    EXPECT(flat(terms), flat(TermList{{0, 0.5}, {1, 1500000000.0}}));
    terms = {{0, 0.5}, {1, 1000000000.0}};
    ConstraintTerms::reduceToIntegers(terms);
    EXPECT(flat(terms), flat(TermList{{0, 1.0}, {1, 2000000000.0}}));
}

void constraintPeriodRules()
{
    // x0 = 0.5 x1: factor 2, period of x1 from 1 to 2 (counted).
    // x2 = 0.25 x1 + x3 + x4 / 3: x1 has period 2, share 0.5, factor 2, stays 2;
    //   x3 is whole; x4 has period 2, share 2/3, factor 3, period 6.
    // A condition without pivot and one whose pivot has period 0 are left out.
    // x3 = 0.3 x0: factor 10, period of x0 from 1 to 10 (counted).
    const std::vector<param::Constraint> conditions{{0, {{1, 0.5}}},
                                                    {2, {{1, 0.25}, {3, 1.0}, {4, 0.3333333333333333}}},
                                                    {-1, {}},
                                                    {5, {{0, 0.5}}},
                                                    {3, {{0, 0.3}}}};
    std::vector<int> periods{1, 1, 1, 1, 2, 0};
    EXPECT(param::ConstraintPeriods::adjust(periods, conditions, 0), 2);
    EXPECT(periods, (std::vector<int>{10, 2, 1, 1, 6, 0}));
    periods = {1, 1, 1, 1, 2, 0};
    EXPECT(param::ConstraintPeriods::adjust(periods, conditions, 4), 1);
    EXPECT(periods, (std::vector<int>{10, 1, 1, 1, 2, 0}));
    periods = {1, 1, 1, 1, 2, 0};
    EXPECT(param::ConstraintPeriods::adjust(periods, conditions, -1), 0);
    EXPECT(periods, (std::vector<int>{1, 1, 1, 1, 2, 0}));

    // No multiple of a share that is not a number or infinite is whole, and
    // 2 * 2147483648.5 is whole but outside the int range: the periods of x1,
    // x2 and x3 stay. x4 has the share 1/3 and the factor 3; the product
    // 3 * 1500000000 is taken in 32 bits (205032704), divided by 3.
    const std::vector<param::Constraint> beyond{{0, {{1, notNumber}, {2, infinite}, {3, 2147483648.5}}},
                                                {0, {{4, 2.2222222222222221e-10}}}};
    periods = {1, 1, 1, 1, 1500000000};
    EXPECT(param::ConstraintPeriods::adjust(periods, beyond, 0), 0);
    EXPECT(periods, (std::vector<int>{1, 1, 1, 1, 68344234}));
}

// Whole test values for the unknowns of a system.
std::vector<double> wholeUnknowns(const param::System& system)
{
    std::vector<double> values(static_cast<std::size_t>(system.unknownCount()));
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = static_cast<double>(static_cast<int>((i * 7 + 3) % 11) - 5);
    }
    return values;
}

// Values of all variables of a system for given values of its unknowns: an
// eliminated variable is the sum its condition names. `open` counts what must
// not occur: a variable that is neither unknown nor pivot of its condition,
// and a term of a condition that is not an unknown.
std::vector<double> variableValues(const param::System& system, const std::vector<double>& unknowns, int& open)
{
    std::vector<double> x(static_cast<std::size_t>(system.variableCount), 0.0);
    open = 0;
    for (int variable = 0; variable < system.variableCount; ++variable) {
        const int unknown = system.variableToUnknown[static_cast<std::size_t>(variable)];
        const int condition = system.constraints.ofVariable(variable);
        if ((unknown == -1) == (condition == -1)) {
            ++open;
        } else if (unknown != -1) {
            x[static_cast<std::size_t>(variable)] = unknowns[static_cast<std::size_t>(unknown)];
        } else if (system.constraints.list()[static_cast<std::size_t>(condition)].pivot != variable) {
            ++open;
        }
    }
    for (int variable = 0; variable < system.variableCount; ++variable) {
        const int condition = system.constraints.ofVariable(variable);
        if (condition == -1) {
            continue;
        }
        double sum = 0.0;
        for (const param::Term& term : system.constraints.list()[static_cast<std::size_t>(condition)].terms) {
            const int unknown = system.variableToUnknown[static_cast<std::size_t>(term.index)];
            if (unknown == -1) {
                ++open;
            } else {
                sum += term.value * unknowns[static_cast<std::size_t>(unknown)];
            }
        }
        x[static_cast<std::size_t>(variable)] = sum;
    }
    return x;
}

// How many coordinates break the transition of a chain for the variable
// values `x`: along a chain, side A minus its first point is side B minus its
// first point turned by the quarter turns of the chain; across a chain with
// offsets, side A is side B turned plus the offsets. Chains that get no
// transition are left out.
int transitionMisfits(const parts::Part& part, const std::vector<double>& x)
{
    static constexpr int cosine[4] = {1, 0, -1, 0};
    static constexpr int sine[4] = {0, 1, 0, -1};
    const param::Chains& chains = *part.chains;
    auto at = [&](int variable) { return x[static_cast<std::size_t>(variable)]; };
    auto differs = [](double a, double b) { return std::fabs(a - b) > 1e-9; };
    int misfits = 0;
    for (int chain = 0; chain < chains.count(); ++chain) {
        const std::size_t c = static_cast<std::size_t>(chain);
        const std::vector<int>& a = chains.sideA[c];
        const std::vector<int>& b = chains.sideB[c];
        if (chains.fullTurnChains.test(c) || chains.curveCheckChains.test(c) || a.empty() || a.size() != b.size()) {
            continue;
        }
        const int quarter = ((static_cast<int>(chains.rotation[c]) % 4) + 4) % 4;
        const double cs = cosine[quarter];
        const double sn = sine[quarter];
        const int offset = chains.offsetVariable[c];
        for (std::size_t k = 0; k < a.size(); ++k) {
            double au = at(2 * a[k]);
            double av = at(2 * a[k] + 1);
            double bu = at(2 * b[k]);
            double bv = at(2 * b[k] + 1);
            if (offset == 0) {
                au -= at(2 * a[0]);
                av -= at(2 * a[0] + 1);
                bu -= at(2 * b[0]);
                bv -= at(2 * b[0] + 1);
            } else {
                au -= at(offset);
                av -= at(offset + 1);
            }
            misfits += differs(au, cs * bu - sn * bv) ? 1 : 0;
            misfits += differs(av, sn * bu + cs * bv) ? 1 : 0;
        }
    }
    return misfits;
}

// A part for the system from a flat mesh of unit squares, each split along
// its rising diagonal into the triangle below it and the one above: frames,
// combed directions and sizes by hand (first direction along x with size
// 0.5, second along y with size 0.25), cut along `cut` (point pairs) with
// `jump` quarter turns left across every cut edge, `singular` as singular
// points and, with `borderCurves`, the border edges as crease curves, as
// well as the edges between the point pairs `creases`. Chains and system
// included.
parts::Part framedPart(parts::Part part, std::initializer_list<std::pair<int, int>> cut, int jump,
                       std::vector<int> singular, bool borderCurves,
                       std::initializer_list<std::pair<int, int>> creases = {})
{
    const Adjacency& a = part.mesh.adjacency();
    const std::size_t faces = static_cast<std::size_t>(part.mesh.faceCount());
    const std::size_t edgeCount = static_cast<std::size_t>(a.edgeCount());

    // Triangle 2i is (p, p+1, p+7) with its reference edge along +x from its
    // first corner; triangle 2i+1 is (p, p+7, p+6) with its reference edge
    // along -y from its third corner.
    field::CrossField cross;
    field::CutGraph cuts;
    for (std::size_t t = 0; t < faces; ++t) {
        const bool lower = t % 2 == 0;
        cross.frameE1.push_back(lower ? Vec3f{1.0f, 0.0f, 0.0f} : Vec3f{0.0f, -1.0f, 0.0f});
        cross.frameE2.push_back(lower ? Vec3f{0.0f, 1.0f, 0.0f} : Vec3f{1.0f, 0.0f, 0.0f});
        cross.startCorner.push_back(static_cast<std::uint8_t>(lower ? 0 : 2));
        cross.referenceLength.push_back(1.0f);
        cuts.combedU.push_back(lower ? std::array<float, 2>{1.0f, 0.0f} : std::array<float, 2>{0.0f, 1.0f});
        cuts.combedV.push_back(lower ? std::array<float, 2>{0.0f, 1.0f} : std::array<float, 2>{-1.0f, 0.0f});
    }
    cross.singularities.points = std::move(singular);
    cross.singularities.index.assign(cross.singularities.points.size(), 1);
    cross.singularities.fullTurn.assign(cross.singularities.points.size(), 0);

    cuts.cutEdgeBits = BitArray(edgeCount);
    cuts.residualJump.assign(edgeCount, 0);
    for (const int e : edgesBetween(a, cut)) {
        cuts.cutEdgeBits.set(static_cast<std::size_t>(e));
        cuts.residualJump[static_cast<std::size_t>(e)] = jump;
    }
    cuts.cutEdgeBits.forEachSet([&](std::size_t e) { cuts.cutEdges.push_back(static_cast<int>(e)); });
    cuts.copy = field::CutCopy::cut(part.mesh, cuts.cutEdges, cuts.cutEdgeBits);
    cuts.curveEdgeKind.assign(edgeCount, field::curve_edge::none);
    if (borderCurves) {
        cuts.curveEdges = a.boundaryEdges();
    }
    for (const int e : edgesBetween(a, creases)) {
        cuts.curveEdges.push_back(e);
    }
    for (const int e : cuts.curveEdges) {
        cuts.curveEdgeKind[static_cast<std::size_t>(e)] = field::curve_edge::crease;
    }
    field::CurveAlignment::orient(part.mesh, cross, cuts);
    part.crossField = std::move(cross);
    part.cutGraph = std::move(cuts);

    field::SizeField sizes;
    sizes.edgeLength = 0.5f;
    sizes.sizeU.assign(faces, 0.5f);
    sizes.sizeV.assign(faces, 0.25f);
    part.sizeField = std::move(sizes);
    part.edgeLength = 0.5f;

    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    param::ChainBuilder::compute(part, param::ChainBuilder::Options{}, progress);
    param::SystemBuilder::compute(part, param::SystemBuilder::Options{}, progress);
    return part;
}

// The flat grid as such a part.
parts::Part flatPart(std::initializer_list<std::pair<int, int>> cut, int jump, std::vector<int> singular,
                     bool borderCurves, std::initializer_list<std::pair<int, int>> creases = {})
{
    parts::Part part;
    flatGrid(part.mesh);
    return framedPart(std::move(part), cut, jump, std::move(singular), borderCurves, creases);
}

// One unit square (points 0 to 3 around it) as such a part: triangle 0 is
// (0, 1, 2), triangle 1 is (0, 2, 3); corner 2 and corner 3 start the
// diagonal they share.
parts::Part squarePart()
{
    parts::Part part;
    std::vector<Vec3f> p{{0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}, {1.0f, 1.0f, 0.0f}, {0.0f, 1.0f, 0.0f}};
    std::vector<int> corners{0, 1, 2, 0, 2, 3};
    part.mesh.swapPoints(p);
    part.mesh.editTopology().swapCorners(corners, FaceMode::Triangles);
    part.mesh.topologyChanged();
    return framedPart(std::move(part), {}, 0, {}, false);
}

// Largest rest |row * x - target| over all rows of a system on the flat grid
// when every copy point gets the coordinates of the perfect grid, u = x / 0.5
// and v = y / 0.25; `misplaced` counts the variables that the conditions put
// elsewhere.
double gridRest(const parts::Part& part, int& misplaced)
{
    const param::System& system = *part.system;
    const std::vector<int>& origin = part.cutGraph->copy.pointOrigins();
    auto gridValue = [&](int variable) {
        const int p = origin[static_cast<std::size_t>(variable / 2)];
        return variable % 2 == 0 ? 2.0 * (p % 6) : 4.0 * (p / 6);
    };
    std::vector<double> unknowns;
    for (const int variable : system.unknownToVariable) {
        unknowns.push_back(gridValue(variable));
    }
    int open = 0;
    const std::vector<double> x = variableValues(system, unknowns, open);
    misplaced = open;
    for (int variable = 0; variable < system.variableCount; ++variable) {
        misplaced += x[static_cast<std::size_t>(variable)] == gridValue(variable) ? 0 : 1;
    }
    double rest = 0.0;
    for (int row = 0; row < system.rows.listCount(); ++row) {
        double sum = -system.targets[static_cast<std::size_t>(row)];
        for (const param::Term& term : system.rows.list(row)) {
            sum += term.value * unknowns[static_cast<std::size_t>(term.index)];
        }
        rest = std::max(rest, std::fabs(sum));
    }
    return rest;
}

// System of an open disc without cut and without singular point. 30 points,
// 40 triangles: 60 variables; the anchor is the first point of the border
// loop (point 0), so 58 unknowns remain and none has to become whole; four
// rows per triangle. The rows of the first triangle (points 0, 1, 7; scales
// 0.5 and 0.25; area factor 1) are, with the anchor terms gone and variables
// 2, 3, 14, 15 as unknowns 0, 1, 12, 13:
//   0.5 u1 = 1 (weight 0.55f),  -0.5 u1 + 0.5 u7 = 0 (weight 1),
//   -0.25 v1 + 0.25 v7 = 1 (weight 0.55f),  0.25 v1 = 0 (weight 1).
// The perfect grid fulfils every row.
// With the border as crease curves the anchor is the corner where curves of
// both directions meet (point 0 again); the four sides become iso lines: the
// lower and the left side run through the anchor (their coordinate is zero),
// the upper and the right side keep one free whole value each. 20 conditions
// (2 anchors, 16 "same as the reference", 2 from the elimination), 40
// unknowns, 2 of them whole.
void systemOnDisc()
{
    using param::TermList;
    {
        const parts::Part part = flatPart({}, 0, {}, false);
        EXPECT(part.system.has_value(), true);
        if (!part.system) {
            return;
        }
        const param::System& system = *part.system;
        EXPECT(system.variableCount, 60);
        EXPECT(system.constraints.count(), 2);
        EXPECT(system.constraints.anchorCount(), 2);
        EXPECT(system.constraints.list()[0].pivot == 0 && system.constraints.list()[1].pivot == 1, true);
        EXPECT(system.unknownCount(), 58);
        EXPECT(system.integerUnknowns.empty() && system.integerPeriods.empty(), true);
        EXPECT(system.energyRowCount, 160);
        EXPECT(system.penaltyRowCount, 0);
        EXPECT(system.rows.listCount(), 160);
        EXPECT(system.rows.otherDimension(), 58);
        EXPECT(system.triangleRowStart.size() == 41 && system.triangleRowStart[1] == 4 &&
                   system.triangleRowStart.back() == 160,
               true);
        EXPECT(flat(system.rows.list(0)), flat(TermList{{0, 0.5}}));
        EXPECT(flat(system.rows.list(1)), flat(TermList{{0, -0.5}, {12, 0.5}}));
        EXPECT(flat(system.rows.list(2)), flat(TermList{{1, -0.25}, {13, 0.25}}));
        EXPECT(flat(system.rows.list(3)), flat(TermList{{1, 0.25}}));
        EXPECT((std::vector<double>(system.targets.begin(), system.targets.begin() + 4)),
               (std::vector<double>{1.0, 0.0, 1.0, 0.0}));
        const double unit = static_cast<double>(1.0f * 0.55f);
        EXPECT((std::vector<double>(system.weights.begin(), system.weights.begin() + 4)),
               (std::vector<double>{unit, 1.0, unit, 1.0}));
        int misplaced = 0;
        EXPECT(gridRest(part, misplaced), 0.0);
        EXPECT(misplaced, 0);
        EXPECT(system.normal.matrix().columns(), 58);
        EXPECT(system.normal.rightHandSide().size(), std::size_t{58});
    }
    {
        const parts::Part part = flatPart({}, 0, {}, true);
        EXPECT(part.system.has_value(), true);
        if (!part.system) {
            return;
        }
        const param::System& system = *part.system;
        EXPECT(part.chains->curves.count(), 4);
        EXPECT(system.constraints.list()[0].pivot == 0 && system.constraints.list()[1].pivot == 1, true);
        // The conditions one by one. The curve chains run along the right
        // side (points 23, 29, 17, 11, 5 in the order of their edges, first
        // coordinate), the left side (18, 24, 12, 6, 0, first coordinate),
        // the upper side (28, 29, 27, 26, 25, 24, second coordinate) and the
        // lower side (4, 5, 3, 2, 1, 0, second coordinate). The first point
        // of each is free and its reference: variables 46, 36, 57 and 9;
        // every other point is tied to it in the order of the points. The
        // anchor ends the left and the lower side: its two rows "reference =
        // 0" are eliminated in the order of their references, 36 before 9,
        // and put into the conditions of these two sides.
        EXPECT(system.curveReference, (std::vector<int>{46, 36, 57, 9}));
        EXPECT(flat(system.constraints.list()),
               flat(std::vector<param::Constraint>{
                   {0, {}}, {1, {}},
                   {58, {{46, 1.0}}}, {34, {{46, 1.0}}}, {22, {{46, 1.0}}}, {10, {{46, 1.0}}},
                   {48, {}}, {24, {}}, {12, {}},
                   {59, {{57, 1.0}}}, {55, {{57, 1.0}}}, {53, {{57, 1.0}}}, {51, {{57, 1.0}}}, {49, {{57, 1.0}}},
                   {11, {}}, {7, {}}, {5, {}}, {3, {}},
                   {36, {}}, {9, {}}}));
        EXPECT(system.sections.anchors, 2);
        EXPECT(system.sections.chainEnds, 2);
        EXPECT(system.sections.isoLines, 18);
        EXPECT(system.constraints.count(), 20);
        EXPECT(system.isoLinesEliminated, true);
        EXPECT(system.penaltyRowCount, 0);
        EXPECT(system.unknownCount(), 40);
        EXPECT(system.integerUnknowns.size(), std::size_t{2});
        int wholeVariables = 0;
        for (const std::uint8_t flag : system.flags) {
            wholeVariables += flag == param::variable_flag::node ? 1 : 0;
        }
        EXPECT(wholeVariables, 22);
        int misplaced = 0;
        EXPECT(gridRest(part, misplaced), 0.0);
        EXPECT(misplaced, 0);
    }
}

// Transitions across a cut chain of four points on the flat grid, for every
// number of quarter turns (also a negative one). Inside the chain the two
// middle points are doubled: 4 conditions.
//  - chain 13-14-15-16 inside the mesh (32 copy points, 64 variables): its
//    ends are not doubled, so without a turn the rows at the ends are empty
//    (rank 0); with a turn they force the last point onto the first
//    (rank 2). The first node stands in as singular point and is the anchor;
//    no unknown has to become whole,
//  - chain 12-13-14-15 from the border to a singular point (33 copy points,
//    66 variables): the border point is doubled, the rows at the ends tie its
//    two copies (rank 2 for every turn); the anchor is the singular point,
//    the second copy of the border node is eliminated, and the two
//    coordinates of its first copy (variables 24, 25) are the unknowns that
//    have to become whole.
// In both, the transition holds at every chain point for whole test values.
// The conditions one by one: in both chains side B of the first middle point
// is copy point 31 (variables 62, 63) and that of the second is copy point 30
// (variables 60, 61); side A is points 15 and 14 for the chain from 16
// (pivots 30, 31, then 28, 29) and points 14 and 13 for the chain from 15
// (pivots 28, 29, then 26, 27). The terms of the first chain point fall away
// (it is the anchor, or its coordinates are eliminated to zero). What is
// left per quarter turn is u_A = cos u_B - sin v_B, v_A = sin u_B + cos v_B.
// At the border the second copy of point 12 (copy point 32, variables 64 and
// 65) lies on side B: its coordinates are the turned ones of the first copy
// (variables 24, 25), u eliminated before v.
void systemTransitions()
{
    using param::Constraint;
    using param::Term;
    // Per quarter turn: u and v of the first middle point on side A, then u
    // and v of the second.
    static const Term middle[4][4] = {{{62, 1.0}, {63, 1.0}, {60, 1.0}, {61, 1.0}},
                                      {{63, -1.0}, {62, 1.0}, {61, -1.0}, {60, 1.0}},
                                      {{62, -1.0}, {63, -1.0}, {60, -1.0}, {61, -1.0}},
                                      {{63, 1.0}, {62, -1.0}, {61, 1.0}, {60, -1.0}}};
    // Per quarter turn: u and v of the second copy of point 12.
    static const Term doubled[4][2] = {{{24, 1.0}, {25, 1.0}},
                                       {{25, 1.0}, {24, -1.0}},
                                       {{24, -1.0}, {25, -1.0}},
                                       {{25, -1.0}, {24, 1.0}}};
    for (const int jump : {0, 1, 2, 3, -1}) {
        const Term* m = middle[(jump + 4) % 4];
        const Term* d = doubled[(jump + 4) % 4];
        {
            const parts::Part part = flatPart({{13, 14}, {14, 15}, {15, 16}}, jump, {}, false);
            EXPECT(part.system.has_value(), true);
            if (!part.system) {
                continue;
            }
            const param::System& system = *part.system;
            const int rank = jump == 0 ? 0 : 2;
            std::vector<Constraint> expected{{26, {}}, {27, {}}, {30, {m[0]}}, {31, {m[1]}}, {28, {m[2]}}, {29, {m[3]}}};
            if (jump != 0) {
                expected.insert(expected.end(), {{32, {}}, {33, {}}});
            }
            EXPECT(flat(system.constraints.list()), flat(expected));
            EXPECT(part.chains->count(), 1);
            EXPECT(system.variableCount, 64);
            EXPECT(system.sections.anchors, 2);
            EXPECT(system.sections.interior, 6);
            EXPECT(system.sections.chainEnds, 6 + rank);
            EXPECT(system.constraints.count(), 6 + rank);
            EXPECT(system.chainEndsEliminated, true);
            EXPECT(system.unknownCount(), 58 - rank);
            EXPECT(system.integerUnknowns.empty(), true);
            int open = 0;
            const std::vector<double> x = variableValues(system, wholeUnknowns(system), open);
            EXPECT(open, 0);
            EXPECT(transitionMisfits(part, x), 0);
        }
        {
            const parts::Part part = flatPart({{12, 13}, {13, 14}, {14, 15}}, jump, {15}, false);
            EXPECT(part.system.has_value(), true);
            if (!part.system) {
                continue;
            }
            const param::System& system = *part.system;
            EXPECT(part.chains->count(), 1);
            EXPECT(system.variableCount, 66);
            EXPECT(system.constraints.list()[0].pivot == 30 && system.constraints.list()[1].pivot == 31, true);
            EXPECT(flat(system.constraints.list()),
                   flat(std::vector<Constraint>{{30, {}}, {31, {}}, {28, {m[0]}}, {29, {m[1]}}, {26, {m[2]}},
                                                {27, {m[3]}}, {64, {d[0]}}, {65, {d[1]}}}));
            EXPECT(system.sections.interior, 6);
            EXPECT(system.sections.chainEnds, 8);
            EXPECT(system.chainEndsEliminated, true);
            EXPECT(system.unknownCount(), 58);
            EXPECT(system.integerUnknowns, (std::vector<int>{24, 25}));
            EXPECT(system.unknownToVariable[24] == 24 && system.unknownToVariable[25] == 25, true);
            int open = 0;
            const std::vector<double> x = variableValues(system, wholeUnknowns(system), open);
            EXPECT(open, 0);
            EXPECT(transitionMisfits(part, x), 0);
        }
    }
}

// The first steps of a system on a part, for checks that look at one step:
// variables and anchors.
struct SystemSteps {
    explicit SystemSteps(parts::Part& part) : context(part, false)
    {
        param::VariableFlags::assign(context, system);
        param::AnchorConstraints::add(context, system);
    }
    param::SystemContext context;
    param::System system;
};

// The transitions of the chain 16-15-14-13 inside the mesh as they are built,
// before the elimination, term by term. Side A is 16, 15, 14, 13, side B is
// 16, 31, 30, 13: the first point (variables 32, 33) is not doubled and free,
// the last one is the anchor (variables 26, 27), whose terms are dropped.
//  - inside, per middle point k and coordinate: the turned coordinate of
//    B[k], minus that of B[0], plus the coordinate of A[0]; sorted by
//    variable, equal variables added up,
//  - at the ends, per coordinate: the turned coordinate of B[3] minus that of
//    B[0], plus the coordinate of A[0], minus that of A[3]. Without a turn
//    everything cancels; the two rows stay in the list without terms.
// Then the cases in which a row names a variable twice, with an odd turn:
// both sides share a point at a place where the cut does not double it. The
// terms are added up there as for an even turn, and what cancels is gone.
// Last a chain point that is the anchor and not doubled: its condition names
// its own pivot, and it is entered for the pivot before the terms of anchors
// are dropped, so the term stays.
void systemTransitionTerms()
{
    using param::Constraint;
    using param::TermList;
    using Rows = std::vector<TermList>;
    struct Turn {
        int jump;
        std::vector<Constraint> interior;
        Rows ends;
    };
    const Turn turns[] = {
        {0, {{30, {{62, 1.0}}}, {31, {{63, 1.0}}}, {28, {{60, 1.0}}}, {29, {{61, 1.0}}}}, {{}, {}}},
        {1,
         {{30, {{32, 1.0}, {33, 1.0}, {63, -1.0}}},
          {31, {{32, -1.0}, {33, 1.0}, {62, 1.0}}},
          {28, {{32, 1.0}, {33, 1.0}, {61, -1.0}}},
          {29, {{32, -1.0}, {33, 1.0}, {60, 1.0}}}},
         {{{32, 1.0}, {33, 1.0}}, {{32, -1.0}, {33, 1.0}}}},
        {2,
         {{30, {{32, 2.0}, {62, -1.0}}},
          {31, {{33, 2.0}, {63, -1.0}}},
          {28, {{32, 2.0}, {60, -1.0}}},
          {29, {{33, 2.0}, {61, -1.0}}}},
         {{{32, 2.0}}, {{33, 2.0}}}},
        {3,
         {{30, {{32, 1.0}, {33, -1.0}, {63, 1.0}}},
          {31, {{32, 1.0}, {33, 1.0}, {62, -1.0}}},
          {28, {{32, 1.0}, {33, -1.0}, {61, 1.0}}},
          {29, {{32, 1.0}, {33, 1.0}, {60, -1.0}}}},
         {{{32, 1.0}, {33, -1.0}}, {{32, 1.0}, {33, 1.0}}}},
    };
    const auto cutPart = [](int jump) { return flatPart({{13, 14}, {14, 15}, {15, 16}}, jump, {}, false); };
    for (const Turn& turn : turns) {
        parts::Part part = cutPart(turn.jump);
        EXPECT(part.chains->sideA[0], (std::vector<int>{16, 15, 14, 13}));
        EXPECT(part.chains->sideB[0], (std::vector<int>{16, 31, 30, 13}));
        SystemSteps steps(part);
        param::TransitionConstraints transitions(steps.context, steps.system);
        transitions.addInterior();
        std::vector<Constraint> expected{{26, {}}, {27, {}}};
        expected.insert(expected.end(), turn.interior.begin(), turn.interior.end());
        EXPECT(flat(steps.system.constraints.list()), flat(expected));
        const Rows ends = transitions.endRows();
        EXPECT(ends.size(), std::size_t{2});
        EXPECT(flat(ends), flat(turn.ends));
    }
    {
        // Side B of the first middle point is the first point of the chain:
        // its turned coordinate minus that of B[0] is gone, the coordinate of
        // A[0] is left. The second middle point is as before.
        parts::Part part = cutPart(1);
        part.chains->sideB[0][1] = 16;
        SystemSteps steps(part);
        param::TransitionConstraints(steps.context, steps.system).addInterior();
        EXPECT(flat(steps.system.constraints.list()),
               flat(std::vector<Constraint>{{26, {}},
                                            {27, {}},
                                            {30, {{32, 1.0}}},
                                            {31, {{33, 1.0}}},
                                            {28, {{32, 1.0}, {33, 1.0}, {61, -1.0}}},
                                            {29, {{32, -1.0}, {33, 1.0}, {60, 1.0}}}}));
    }
    {
        // Side B ends at the point it starts from: the turned coordinates
        // cancel, the coordinate of A[0] is left (that of A[3] is an anchor).
        parts::Part part = cutPart(1);
        part.chains->sideB[0].back() = 16;
        SystemSteps steps(part);
        EXPECT(flat(param::TransitionConstraints(steps.context, steps.system).endRows()),
               flat(Rows{{{32, 1.0}}, {{33, 1.0}}}));
    }
    {
        // Side A ends at the point it starts from: its two coordinates
        // cancel, the turned ones of side B are left: -v of B[3] (an anchor)
        // + v of B[0] for the first row, u of B[3] (an anchor) - u of B[0]
        // for the second.
        parts::Part part = cutPart(1);
        part.chains->sideA[0].back() = 16;
        SystemSteps steps(part);
        EXPECT(flat(param::TransitionConstraints(steps.context, steps.system).endRows()),
               flat(Rows{{{33, 1.0}}, {{32, -1.0}}}));
    }
    {
        // Both sides closed: nothing is left of either row.
        parts::Part part = cutPart(1);
        part.chains->sideA[0].back() = 16;
        part.chains->sideB[0].back() = 16;
        SystemSteps steps(part);
        EXPECT(flat(param::TransitionConstraints(steps.context, steps.system).endRows()), flat(Rows{{}, {}}));
    }
    {
        // The first middle point is the anchor on both sides, no turn: u of
        // B[1] is the pivot itself, the coordinates of B[0] and A[0] cancel.
        // The pivot stays as a term, the entries of the anchor point to the
        // new conditions, and the two anchor conditions stay in the list.
        parts::Part part = cutPart(0);
        part.chains->sideA[0][1] = 13;
        part.chains->sideB[0][1] = 13;
        SystemSteps steps(part);
        param::TransitionConstraints(steps.context, steps.system).addInterior();
        const param::ConstraintSet& set = steps.system.constraints;
        EXPECT(flat(set.list()), flat(std::vector<Constraint>{{26, {}},
                                                              {27, {}},
                                                              {26, {{26, 1.0}}},
                                                              {27, {{27, 1.0}}},
                                                              {28, {{60, 1.0}}},
                                                              {29, {{61, 1.0}}}}));
        EXPECT(set.ofVariable(26) == 2 && set.ofVariable(27) == 3 && set.anchorCount() == 2, true);
    }
    // A turn that is no number or too large for a whole number is no turn; a
    // negative one counts backwards.
    for (const float turn : {std::numeric_limits<float>::quiet_NaN(), 3e9f, -3e9f}) {
        const param::TransitionConstraints::Turn none = param::TransitionConstraints::turnOf(turn);
        EXPECT(!none.odd && none.cosine == 1.0 && none.sine == 0.0, true);
    }
    const param::TransitionConstraints::Turn back = param::TransitionConstraints::turnOf(-1.0f);
    EXPECT(back.odd && back.cosine == 0.0 && back.sine == -1.0, true);
}

// Entry of a matrix in compressed columns, 0 where none is stored.
double entryOf(const solver::CompactMatrix& matrix, int row, int column)
{
    for (int k = matrix.columnStart(column); k < matrix.columnEnd(column); ++k) {
        if (matrix.row(k) == row) {
            return matrix.value(k);
        }
    }
    return 0.0;
}

// An elimination of iso line rows that fails, on the open disc without cut
// (60 variables, all continuous, anchor = variables 0 and 1). The rows are
//   2 x40 + x20 = 0,   t x30 + 2 x50 = 0 (t = 2^-34),   2 x44 + 4 x24 = 0,
//   4 x44 + x10 = 0.
// Columns in order of occurrence: x40, x20, x30, x50, x44, x24, x10.
//  - column 0 has no leading 1; the 1 of x20 is exchanged into it and row 0,
//    x20 + 2 x40, is the first pivot row,
//  - column 1 (now x40): the 1 of x10 in the last row is exchanged into it,
//    and that row, x10 + 4 x44, becomes row 1; x40 now has the last column,
//  - column 2 (x30) has only t: the row t x30 + 2 x50 becomes row 2, the
//    first pivot that is not 1,
//  - column 3 (x50) leads in no row; column 4 (x44): 2 x44 + 4 x24.
// From the last pivot backwards: x44 = -2 x24; then t is too small to divide
// by and the run fails. The rows stay as they are, and they and the reached
// condition speak of variables again; the three conditions that were not
// reached are appended empty and are the condition of no variable.
// The rows that are left are rows of the system behind the 160 energy rows,
// with the conditions put in and weight 1e20: x44 is replaced by -2 x24, so
// the last row is empty. Variables 0, 1 and 44 are no unknowns: variable v is
// unknown v - 2 below 44 and v - 3 above. Every penalty row adds weight times
// the product of its coefficients to the normal matrix (exactly: all factors
// are powers of two times 1e20) and nothing to the right-hand side, its
// target being 0. The same set-up without the penalty rows gives the rest.
void systemPenaltyRows()
{
    using param::Constraint;
    using param::TermList;
    using Rows = std::vector<TermList>;
    const double tiny = 0x1p-34;
    const Rows isoRows{{{40, 2.0}, {20, 1.0}}, {{30, tiny}, {50, 2.0}}, {{44, 2.0}, {24, 4.0}}, {{44, 4.0}, {10, 1.0}}};
    QuietProgress q;

    parts::Part part = flatPart({}, 0, {}, false);
    param::SystemSetup setup(part, false, q.progress);
    setup.constrainChains();
    Rows rows = isoRows;
    EXPECT(setup.eliminate(rows, param::EliminationMode::IsoLines, {}, false), false);
    EXPECT(flat(rows),
           flat(Rows{{{20, 1.0}, {40, 2.0}}, {{10, 1.0}, {44, 4.0}}, {{30, tiny}, {50, 2.0}}, {{44, 2.0}, {24, 4.0}}}));
    const param::ConstraintSet& set = setup.system().constraints;
    EXPECT(flat(set.list()),
           flat(std::vector<Constraint>{{0, {}}, {1, {}}, {-1, {}}, {-1, {}}, {-1, {}}, {44, {{24, -2.0}}}}));
    std::vector<int> pivots;
    for (int variable = 0; variable < setup.system().variableCount; ++variable) {
        if (set.ofVariable(variable) != -1) {
            pivots.insert(pivots.end(), {variable, set.ofVariable(variable)});
        }
    }
    EXPECT(pivots, (std::vector<int>{0, 0, 1, 1, 44, 5}));

    // The same conditions without penalty rows.
    parts::Part plainPart = flatPart({}, 0, {}, false);
    param::SystemSetup plain(plainPart, false, q.progress);
    plain.constrainChains();
    Rows unused = isoRows;
    EXPECT(plain.eliminate(unused, param::EliminationMode::IsoLines, {}, false), false);
    EXPECT(plain.assemble({}) == param::SystemSetup::Outcome::Built, true);
    EXPECT(setup.assemble(rows) == param::SystemSetup::Outcome::Built, true);
    const param::System& system = setup.system();
    const param::System& bare = plain.system();

    EXPECT(system.unknownCount(), 57);
    EXPECT(system.energyRowCount, 160);
    EXPECT(system.penaltyRowCount, 4);
    EXPECT(system.rows.listCount(), 164);
    EXPECT(bare.rows.listCount(), 160);
    EXPECT(flat(system.rows.list(160)), flat(TermList{{18, 1.0}, {38, 2.0}}));
    EXPECT(flat(system.rows.list(161)), flat(TermList{{8, 1.0}, {22, -8.0}}));
    EXPECT(flat(system.rows.list(162)), flat(TermList{{28, tiny}, {47, 2.0}}));
    EXPECT(system.rows.list(163).empty(), true);
    EXPECT((std::vector<double>(system.targets.begin() + 160, system.targets.end())),
           (std::vector<double>(4, 0.0)));
    EXPECT((std::vector<double>(system.weights.begin() + 160, system.weights.end())),
           (std::vector<double>(4, 1e20)));

    // Normal matrix: unknowns (i, j) of a penalty row and weight times the
    // product of their coefficients.
    struct Product {
        int i;
        int j;
        double value;
    };
    const Product products[] = {{18, 18, 1e20},  {18, 38, 2e20},          {38, 38, 4e20},
                                {8, 8, 1e20},    {8, 22, -8e20},          {22, 22, 64e20},
                                {28, 28, 1e20 * tiny * tiny}, {28, 47, 2e20 * tiny}, {47, 47, 4e20}};
    const solver::CompactMatrix& matrix = system.normal.matrix();
    const solver::CompactMatrix& bareMatrix = bare.normal.matrix();
    for (const Product& p : products) {
        EXPECT(entryOf(matrix, p.i, p.j), entryOf(bareMatrix, p.i, p.j) + p.value);
        EXPECT(entryOf(matrix, p.j, p.i), entryOf(bareMatrix, p.j, p.i) + p.value);
    }
    // The large products swallow what the energy rows add at their places.
    EXPECT(entryOf(matrix, 18, 18), 1e20);
    EXPECT(entryOf(matrix, 38, 18), 2e20);
    EXPECT(entryOf(matrix, 38, 38), 4e20);
    EXPECT(entryOf(matrix, 22, 8), -8e20);
    EXPECT(entryOf(matrix, 22, 22), 64e20);
    EXPECT(entryOf(matrix, 47, 28), 2e20 * tiny);
    // Everywhere else the matrix is the one without penalty rows: 6 places on
    // the diagonal and 3 pairs differ; the pairs are points that share no
    // triangle, so their 6 entries are new.
    int different = 0;
    for (int column = 0; column < matrix.columns(); ++column) {
        for (int k = matrix.columnStart(column); k < matrix.columnEnd(column); ++k) {
            different += matrix.value(k) != entryOf(bareMatrix, matrix.row(k), column) ? 1 : 0;
        }
    }
    EXPECT(different, 12);
    EXPECT(matrix.nonZeros(), bareMatrix.nonZeros() + 6);
    EXPECT(system.normal.rightHandSide() == bare.normal.rightHandSide(), true);
    EXPECT(system.normal.rightHandSide().size(), std::size_t{57});
}

// A row term without an unknown: with triangle 0 marked as mirrored half and
// symmetry switched off, the points 0, 1 and 7 (whose first corner lies in
// that triangle) get no unknowns, but every triangle keeps its energy rows.
// The rows of the triangles around these points name their variables: there
// is no system, and the answer says why.
void systemWithoutUnknown()
{
    parts::Part part = flatPart({}, 0, {}, false);
    EXPECT(part.system.has_value(), true);
    part.cutGraph->mirrorHalf = BitArray(static_cast<std::size_t>(part.mesh.faceCount()));
    part.cutGraph->mirrorHalf.set(0);
    QuietProgress q;
    EXPECT(param::SystemBuilder::compute(part, param::SystemBuilder::Options{}, q.progress) ==
               param::SystemSetup::Outcome::TermWithoutUnknown,
           true);
    EXPECT(part.system.has_value(), false);
}

// System of a closed part with cross field and cut graph from the core: every
// variable is an unknown or the pivot of its condition, conditions name
// unknowns only, both eliminations succeed, and the transitions of all chains
// hold for whole test values. A torus whose cuts turn the coordinates by
// whole turns only has two coordinates per point of the mesh, one anchor and
// two periods: twice its points plus 2 unknowns.
void systemOnClosedPart(parts::Part& part, bool torus)
{
    if (!part.chains) {
        return;
    }
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    param::SystemBuilder::compute(part, param::SystemBuilder::Options{}, progress);
    EXPECT(part.system.has_value(), true);
    if (!part.system) {
        return;
    }
    const param::System& system = *part.system;
    EXPECT(system.variableCount,
           2 * part.cutGraph->copy.pointCount() + 2 * static_cast<int>(part.chains->genusChains.size()));
    EXPECT(system.constraints.anchorCount(), 2);
    EXPECT(system.chainEndsEliminated && system.isoLinesEliminated, true);
    EXPECT(system.unresolvedTerms, 0);
    EXPECT(system.rows.listCount(), 4 * part.mesh.faceCount());
    int open = 0;
    const std::vector<double> x = variableValues(system, wholeUnknowns(system), open);
    EXPECT(open, 0);
    EXPECT(transitionMisfits(part, x), 0);
    bool ascending = true;
    for (std::size_t i = 1; i < system.integerUnknowns.size(); ++i) {
        ascending = ascending && system.integerUnknowns[i - 1] < system.integerUnknowns[i];
    }
    EXPECT(ascending, true);
    if (torus) {
        bool turned = false;
        for (const float turns : part.chains->rotation) {
            turned = turned || static_cast<int>(turns) % 4 != 0;
        }
        EXPECT(turned, false);
        EXPECT(system.unknownCount(), 2 * part.mesh.pointCount() + 2);
    }
}

// A rounding step that finds no values.
class NoValues final : public param::RoundingStep {
public:
    Outcome select(const param::SolveState&, solver::DriverRequest&) override
    {
        ++calls;
        return Outcome::Failed;
    }
    int calls = 0;
};

// Solves the system of a part without the rounding by regions: a step that
// finds no values stands in for it, so the driver rounds one unknown per
// round.
void solvePart(parts::Part& part)
{
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    NoValues regions;
    param::BorderRounding borders;
    param::SolutionBuilder::compute(part, param::SolutionBuilder::Options{}, progress,
                                    param::SolutionBuilder::Steps{regions, borders});
}

// Largest distance of the coordinates of a solved flat grid from the perfect
// grid, u = x / 0.5 and v = y / 0.25 at every copy point.
double gridMisfit(const parts::Part& part)
{
    const std::vector<int>& origin = part.cutGraph->copy.pointOrigins();
    double worst = 0.0;
    for (std::size_t w = 0; w < part.solution->uv.size(); ++w) {
        const int p = origin[w];
        worst = std::max(worst, std::fabs(static_cast<double>(part.solution->uv[w][0]) - 2.0 * (p % 6)));
        worst = std::max(worst, std::fabs(static_cast<double>(part.solution->uv[w][1]) - 4.0 * (p / 6)));
    }
    return worst;
}

// How many coordinates of a solved part that had to become whole numbers are
// not exactly whole.
int notWhole(const parts::Part& part)
{
    const param::System& system = *part.system;
    int count = 0;
    for (int variable = 0; variable < 2 * system.copyPointCount; ++variable) {
        if (system.flags[static_cast<std::size_t>(variable)] == param::variable_flag::continuous) {
            continue;
        }
        const float value = part.solution->uv[static_cast<std::size_t>(variable / 2)][static_cast<std::size_t>(variable % 2)];
        count += value == std::round(value) ? 0 : 1;
    }
    return count;
}

// Largest break of a transition in the coordinates of a solved part: along a
// chain, side A minus its first point is side B minus its first point turned
// by the quarter turns of the chain. Chains that get no transition are left out.
double chainBreak(const parts::Part& part)
{
    static constexpr int cosine[4] = {1, 0, -1, 0};
    static constexpr int sine[4] = {0, 1, 0, -1};
    const param::Chains& chains = *part.chains;
    const param::UvOutput::Coordinates& uv = part.solution->uv;
    auto at = [&](int w, std::size_t axis) { return static_cast<double>(uv[static_cast<std::size_t>(w)][axis]); };
    double worst = 0.0;
    for (int chain = 0; chain < chains.count(); ++chain) {
        const std::size_t c = static_cast<std::size_t>(chain);
        const std::vector<int>& a = chains.sideA[c];
        const std::vector<int>& b = chains.sideB[c];
        if (chains.fullTurnChains.test(c) || chains.curveCheckChains.test(c) || a.empty() || a.size() != b.size()) {
            continue;
        }
        const int quarter = ((static_cast<int>(chains.rotation[c]) % 4) + 4) % 4;
        const double cs = cosine[quarter];
        const double sn = sine[quarter];
        for (std::size_t k = 1; k < a.size(); ++k) {
            const double au = at(a[k], 0) - at(a[0], 0);
            const double av = at(a[k], 1) - at(a[0], 1);
            const double bu = at(b[k], 0) - at(b[0], 0);
            const double bv = at(b[k], 1) - at(b[0], 1);
            worst = std::max(worst, std::fabs(au - (cs * bu - sn * bv)));
            worst = std::max(worst, std::fabs(av - (sn * bu + cs * bv)));
        }
    }
    return worst;
}

// Solution on the flat grid: frames along x and y with sizes 0.5 and 0.25 ask
// for u = 2x and v = 4y, which fulfils every row.
//  - without curves no unknown has to become whole: the system is solved
//    once and no round is asked,
//  - with the border as crease curves two unknowns have to become whole (u
//    of the right side, v of the upper side): one round of reweighting, then
//    one unknown per round because the regions find no values; the corner
//    opposite the anchor lands exactly on (10, 16).
void solutionOnDisc()
{
    {
        parts::Part part = flatPart({}, 0, {}, false);
        solvePart(part);
        EXPECT(part.solution.has_value(), true);
        if (!part.solution) {
            return;
        }
        EXPECT(part.solution->uv.size(), std::size_t{30});
        EXPECT(part.solution->rounds, 0);
        EXPECT(part.solution->integerUnknowns, 0);
        EXPECT(gridMisfit(part) < 1e-4, true);
    }
    {
        parts::Part part = flatPart({}, 0, {}, true);
        solvePart(part);
        EXPECT(part.solution.has_value(), true);
        if (!part.solution) {
            return;
        }
        EXPECT(part.solution->rounds, 3);
        EXPECT(part.solution->integerUnknowns, 2);
        EXPECT(gridMisfit(part) < 1e-4, true);
        EXPECT(notWhole(part), 0);
        EXPECT(part.solution->uv[29] == (std::array<float, 2>{10.0f, 16.0f}), true);
    }
}

// Stands in for the rounding by regions: on its first call it delivers the
// first position of the integer list, rounded.
class FirstPosition final : public param::RoundingStep {
public:
    Outcome select(const param::SolveState& state, solver::DriverRequest& request) override
    {
        ++calls;
        const int unknown = state.driver.integerUnknowns().front();
        request.selection.push_back({0, std::round((*request.x)[static_cast<std::size_t>(unknown)])});
        return Outcome::Delivered;
    }
    int calls = 0;
};

// Order of the rounding steps on the flat grid with border curves (two
// unknowns to round): reweighting, then the regions deliver one position,
// then everything left is rounded at once - three rounds, the regions asked
// once.
void solutionRounds()
{
    parts::Part part = flatPart({}, 0, {}, true);
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    FirstPosition regions;
    FirstPosition borders;
    param::SolutionBuilder::compute(part, param::SolutionBuilder::Options{}, progress,
                                    param::SolutionBuilder::Steps{regions, borders});
    EXPECT(part.solution.has_value(), true);
    if (!part.solution) {
        return;
    }
    EXPECT(regions.calls, 1);
    EXPECT(borders.calls, 0);
    EXPECT(part.solution->rounds, 3);
    EXPECT(gridMisfit(part) < 1e-4, true);
    EXPECT(notWhole(part), 0);
}

// Solution across a cut chain from the border to a singular point, for every
// number of quarter turns: the copies of the chain ends get whole
// coordinates and the transition holds at every chain point.
void solutionTransitions()
{
    for (const int jump : {0, 1, 2, 3, -1}) {
        parts::Part part = flatPart({{12, 13}, {13, 14}, {14, 15}}, jump, {15}, false);
        solvePart(part);
        EXPECT(part.solution.has_value(), true);
        if (!part.solution) {
            continue;
        }
        EXPECT(part.solution->integerUnknowns, 2);
        EXPECT(part.solution->rounds, 3);
        EXPECT(notWhole(part), 0);
        EXPECT(chainBreak(part) < 1e-3, true);
        EXPECT(part.solution->singularities.points, (std::vector<int>{15}));
    }
}

// Solution of a closed part: singular points and nodes of the cuts lie on
// whole coordinates, the transitions of all chains hold, and without the
// regions every unknown to round takes one round after the reweighting.
void solutionOnClosedPart(parts::Part& part)
{
    if (!part.system) {
        return;
    }
    const int toRound = static_cast<int>(part.system->integerUnknowns.size());
    solvePart(part);
    EXPECT(part.solution.has_value(), true);
    if (!part.solution) {
        return;
    }
    EXPECT(toRound > 0, true);
    EXPECT(part.solution->integerUnknowns, toRound);
    EXPECT(part.solution->rounds, toRound + 1);
    EXPECT(notWhole(part), 0);
    EXPECT(chainBreak(part) < 1e-3, true);
}

// Values of the unknowns of a flat grid system when copy point (x, y) gets
// u = ux * x + uy * y and v = 4 y.
std::vector<double> shearedGrid(const parts::Part& part, double ux, double uy)
{
    const std::vector<int>& origin = part.cutGraph->copy.pointOrigins();
    std::vector<double> unknowns;
    for (const int variable : part.system->unknownToVariable) {
        const int p = origin[static_cast<std::size_t>(variable / 2)];
        unknowns.push_back(variable % 2 == 0 ? ux * (p % 6) + uy * (p / 6) : 4.0 * (p / 6));
    }
    return unknowns;
}

// Largest relative distance of the weights of rows `first`, `first + 4`, ...
// from `factor` times their old value.
double weightMisfit(const std::vector<double>& now, const std::vector<double>& before, std::size_t first, double factor)
{
    double worst = 0.0;
    for (std::size_t i = first; i < now.size(); i += 4) {
        worst = std::max(worst, std::fabs(now[i] / before[i] - factor));
    }
    return worst;
}

// Reweighting on the flat grid (sizes 0.5 and 0.25, every triangle alike, so
// averaging keeps the factors):
//  - u stretched to twice the step (u = x): the step 1 lies between 1.8 and
//    4 times the size, the first row of every triangle gets
//    1 + 3 * (1 - 0.9) / (2 - 0.9); the other rows stay,
//  - u turned over (u = -2x): all four rows get 1 + 5,
//  - u sheared (u = 2x + 2y): the grid directions are 45 degrees off a right
//    angle, the two rows across get 1 + (3 * (45 - 15) / 75 + 1).
void reweightingFactors()
{
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    struct Case {
        double ux;
        double uy;
        std::array<double, 4> factors;
    };
    const Case cases[] = {{1.0, 0.0, {1.0 + 3.0 * (1.0 - 0.9) / (2.0 - 0.9), 1.0, 1.0, 1.0}},
                          {-2.0, 0.0, {6.0, 6.0, 6.0, 6.0}},
                          {2.0, 2.0, {1.0, 3.2, 1.0, 3.2}}};
    for (const Case& c : cases) {
        parts::Part part = flatPart({}, 0, {}, false);
        EXPECT(part.system.has_value(), true);
        if (!part.system) {
            continue;
        }
        const std::vector<double> before = part.system->weights;
        const param::SystemContext context(part, false);
        param::Reweighting::apply(context, *part.system, shearedGrid(part, c.ux, c.uy), progress);
        for (std::size_t row = 0; row < 4; ++row) {
            EXPECT(weightMisfit(part.system->weights, before, row, c.factors[row]) < 1e-5, true);
        }
    }
}

// Mirror corners of a triangle: the mirror triangle runs the other way
// round, whichever of its corners comes first; a triangle that runs the same
// way or has no mirror triangle gets none.
void mirrorCorners()
{
    const std::vector<int> pointMirror{3, 4, 5, 0, 1, 2};
    const std::array<int, 3> none{-1, -1, -1};
    auto corners = [&](std::initializer_list<int> second, std::vector<int> faceMirror) {
        Mesh m;
        std::vector<Vec3f> p(6);
        std::vector<int> c{0, 1, 2};
        c.insert(c.end(), second);
        m.swapPoints(p);
        m.editTopology().swapCorners(c, FaceMode::Triangles);
        m.topologyChanged();
        return param::UvOutput::mirroredCorners(m.topology(), pointMirror, faceMirror, 1);
    };
    EXPECT(corners({3, 5, 4}, {1, 0}) == (std::array<int, 3>{0, 2, 1}), true);
    EXPECT(corners({4, 3, 5}, {1, 0}) == (std::array<int, 3>{1, 0, 2}), true);
    EXPECT(corners({5, 4, 3}, {1, 0}) == (std::array<int, 3>{2, 1, 0}), true);
    EXPECT(corners({3, 4, 5}, {1, 0}) == none, true);
    EXPECT(corners({3, 5, 4}, {1, -1}) == none, true);
}

// Start points of the extraction on the flat grid.
//  - Without singular point every border point with whole coordinates
//    becomes a seed, once: with the perfect grid all 18 border points, with
//    one coordinate off 17; a part that already has a seed gets none unless
//    the borders are to be frozen.
//  - Symmetry (mirror x -> 5 - x, the triangles right of x = 3 mirrored): a
//    point that touches only mirrored triangles is dropped, every other one
//    is followed by its mirror point.
void solutionSeeds()
{
    parts::Part part = flatPart({}, 0, {}, false);
    const param::SystemContext context(part, false);
    param::Solution perfect;
    for (int p = 0; p < 30; ++p) {
        perfect.uv.push_back({2.0f * static_cast<float>(p % 6), 4.0f * static_cast<float>(p / 6)});
    }
    {
        param::Solution s = perfect;
        param::SolutionSeeds::addBorderPoints(context, false, s);
        EXPECT(s.seeds.size(), std::size_t{18});
        std::vector<int> sorted = s.seeds;
        std::sort(sorted.begin(), sorted.end());
        EXPECT(std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end(), true);
        EXPECT(std::binary_search(sorted.begin(), sorted.end(), 7), false);
    }
    {
        param::Solution s = perfect;
        s.uv[1][0] = 2.5f;
        param::SolutionSeeds::addBorderPoints(context, false, s);
        EXPECT(s.seeds.size(), std::size_t{17});
        EXPECT(std::find(s.seeds.begin(), s.seeds.end(), 1) == s.seeds.end(), true);
    }
    {
        param::Solution s = perfect;
        s.seeds = {7};
        param::SolutionSeeds::addBorderPoints(context, false, s);
        EXPECT(s.seeds, (std::vector<int>{7}));
        param::SolutionSeeds::addBorderPoints(context, true, s);
        EXPECT(s.seeds.size(), std::size_t{19});
    }
    {
        std::vector<int> pointMirror(30);
        for (int p = 0; p < 30; ++p) {
            pointMirror[static_cast<std::size_t>(p)] = 6 * (p / 6) + (5 - p % 6);
        }
        const prep::SymmetryCache pairs = prep::SymmetryCache::fromPoints(pointMirror, BitArray(30), std::nullopt);
        field::CutGraph& graph = *part.cutGraph;
        graph.mirrorHalf = BitArray(40);
        for (int j = 0; j < 4; ++j) {
            for (int i = 3; i < 5; ++i) {
                graph.mirrorHalf.set(static_cast<std::size_t>(2 * (5 * j + i)));
                graph.mirrorHalf.set(static_cast<std::size_t>(2 * (5 * j + i) + 1));
            }
        }
        graph.mirrorHalfCount = 16;
        param::Solution s;
        s.singularities.points = {7, 10, 9};
        s.singularities.index = {1, 3, 2};
        s.singularities.fullTurn = {0, 0, 1};
        s.seeds = {28, 0};
        param::SolutionSeeds::mirror(context, pairs, s);
        EXPECT(s.singularities.points, (std::vector<int>{7, 10, 9, 8}));
        EXPECT(s.singularities.index, (std::vector<int>{1, 1, 2, 2}));
        EXPECT(s.singularities.fullTurn, (std::vector<std::uint8_t>{0, 0, 1, 1}));
        EXPECT(s.seeds, (std::vector<int>{0, 5}));
        EXPECT(s.singularities.indexSum, 6);
    }
}

// From which neighbour the averaging of the reweighting and the smoothing of
// the sizes take over, and what. Mesh: the unit square (triangle 0 with the
// corners 0 to 2, triangle 1 with the corners 3 to 5) and a third triangle
// (3, 3, 0) without area.
//  - corners 2 and 3 start the diagonal: each has the other triangle as
//    neighbour,
//  - a border edge gives no neighbour; corner 6 runs from point 3 to point 3
//    and has no edge at all; the edge 0-3 joins three corners (5, 7 and 8)
//    and gives none of them a neighbour,
//  - across an odd number of quarter turns, also a negative one, the two
//    directions change places,
//  - a crease limits in both uses, a curve edge of an edge loop only where
//    every curve limits; the bit of the corner picks the direction that
//    passes, and a corner beyond the bits lets the second one pass.
void smoothingNeighbourRules()
{
    Mesh m;
    std::vector<Vec3f> p{{0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}, {1.0f, 1.0f, 0.0f}, {0.0f, 1.0f, 0.0f}};
    std::vector<int> corners{0, 1, 2, 0, 2, 3, 3, 3, 0};
    m.swapPoints(p);
    m.editTopology().swapCorners(corners, FaceMode::Triangles);
    m.topologyChanged();
    const Adjacency& a = m.adjacency();
    EXPECT(a.cornerEdge(6), -1);
    EXPECT(a.edgeCorners(a.cornerEdge(5)).size(), std::size_t{3});
    const auto diagonal = static_cast<std::size_t>(a.edgeBetween(0, 2));

    field::CutGraph cuts;
    cuts.curveEdgeKind.assign(static_cast<std::size_t>(a.edgeCount()), field::curve_edge::none);
    cuts.residualJump.assign(static_cast<std::size_t>(a.edgeCount()), 0);
    cuts.curveFollowsU = BitArray(3);
    cuts.curveFollowsU.set(2);
    using Rule = field::SmoothingNeighbours;
    const Rule every(a, cuts, Rule::Limit::EveryCurve);
    const Rule creases(a, cuts, Rule::Limit::Creases);
    // Neighbour triangle, directions change places, first passes, second
    // passes; all -1 without neighbour.
    using Answer = std::vector<int>;
    auto answer = [](const Rule& rule, int corner) {
        const auto n = rule.across(corner);
        return n ? Answer{n->triangle, n->odd ? 1 : 0, n->take.u ? 1 : 0, n->take.v ? 1 : 0} : Answer{-1, -1, -1, -1};
    };
    for (const Rule* rule : {&every, &creases}) {
        EXPECT(answer(*rule, 2), (Answer{1, 0, 1, 1}));
        EXPECT(answer(*rule, 3), (Answer{0, 0, 1, 1}));
        for (const int corner : {0, 1, 4, 5, 6, 7, 8}) {
            EXPECT(answer(*rule, corner), (Answer{-1, -1, -1, -1}));
        }
    }
    for (const int jump : {1, 3, -1}) {
        cuts.residualJump[diagonal] = jump;
        EXPECT(answer(every, 2), (Answer{1, 1, 1, 1}));
        EXPECT(answer(creases, 3), (Answer{0, 1, 1, 1}));
    }
    cuts.residualJump[diagonal] = 2;
    EXPECT(answer(every, 2), (Answer{1, 0, 1, 1}));

    cuts.curveEdgeKind[diagonal] = field::curve_edge::edgeLoop;
    EXPECT(answer(every, 2), (Answer{1, 0, 1, 0}));
    EXPECT(answer(every, 3), (Answer{0, 0, 0, 1}));
    EXPECT(answer(creases, 2), (Answer{1, 0, 1, 1}));
    EXPECT(answer(creases, 3), (Answer{0, 0, 1, 1}));
    cuts.curveEdgeKind[diagonal] = field::curve_edge::crease;
    for (const Rule* rule : {&every, &creases}) {
        EXPECT(answer(*rule, 2), (Answer{1, 0, 1, 0}));
        EXPECT(answer(*rule, 3), (Answer{0, 0, 0, 1}));
    }
}

// Reweighting with its averaging over neighbours, on one unit square with
// u = x and v = 4y. Triangle 0 asks for the sizes 1 and 0.25 and is mapped
// well: factors 1, 1, 1, 1. Triangle 1 asks for 0.125 along its first
// direction; the step 1 is more than 4 times that, and its first row gets
// 1 + 3: factors 4, 1, 1, 1. The two triangles are neighbours across the
// diagonal and have no other neighbour. After the ten rounds of averaging:
//  - plain diagonal: both hold the mean of the two, 2.5 on the first row,
//  - an odd number of quarter turns left across the diagonal: the rows of
//    the two directions change places on the way. The first row of triangle
//    1 meets the third of triangle 0: 2.5 on these two, and the first row of
//    triangle 0 stays 1,
//  - the diagonal as a crease that follows the first direction at the corner
//    of triangle 0 and the second at the corner of triangle 1: triangle 0
//    takes over the rows of the first direction only, triangle 1 those of
//    the second only. The first row of triangle 1 stays 4, and that of
//    triangle 0 comes half the way closer to it in every round: 4 - 3 / 2^10,
//  - the diagonal as curve edge of an edge loop limits nothing here: as the
//    plain diagonal.
void reweightingAveraging()
{
    QuietProgress q;
    struct Case {
        int jump;
        std::uint8_t curve;
        std::array<double, 4> first;
        std::array<double, 4> second;
    };
    const double nearly4 = 4.0 - 3.0 / 1024.0;
    const Case cases[] = {
        {0, field::curve_edge::none, {2.5, 1.0, 1.0, 1.0}, {2.5, 1.0, 1.0, 1.0}},
        {1, field::curve_edge::none, {1.0, 1.0, 2.5, 1.0}, {2.5, 1.0, 1.0, 1.0}},
        {0, field::curve_edge::crease, {nearly4, 1.0, 1.0, 1.0}, {4.0, 1.0, 1.0, 1.0}},
        {0, field::curve_edge::edgeLoop, {2.5, 1.0, 1.0, 1.0}, {2.5, 1.0, 1.0, 1.0}},
    };
    for (const Case& c : cases) {
        parts::Part part = squarePart();
        EXPECT(part.system.has_value(), true);
        if (!part.system) {
            continue;
        }
        param::System& system = *part.system;
        EXPECT(system.triangleRowStart, (std::vector<int>{0, 4, 8}));
        part.sizeField->sizeU = {1.0f, 0.125f};
        field::CutGraph& graph = *part.cutGraph;
        const auto diagonal = static_cast<std::size_t>(part.mesh.adjacency().edgeBetween(0, 2));
        graph.residualJump[diagonal] = c.jump;
        graph.curveEdgeKind[diagonal] = c.curve;
        graph.curveFollowsU = BitArray(6);
        graph.curveFollowsU.set(2);

        const std::vector<int>& origin = graph.copy.pointOrigins();
        std::vector<double> x;
        for (const int variable : system.unknownToVariable) {
            const Vec3f& point = part.mesh.points()[static_cast<std::size_t>(origin[static_cast<std::size_t>(variable / 2)])];
            x.push_back(variable % 2 == 0 ? static_cast<double>(point.x) : 4.0 * static_cast<double>(point.y));
        }
        const std::vector<double> before = system.weights;
        const param::SystemContext context(part, false);
        param::Reweighting::apply(context, system, x, q.progress);
        std::vector<double> expected;
        for (std::size_t row = 0; row < 8; ++row) {
            expected.push_back(before[row] * (row < 4 ? c.first[row] : c.second[row - 4]));
        }
        EXPECT(system.weights, expected);
    }
}

// A stop request that arrives in a round ends the run before the round does
// anything: no step is asked, nothing is reported, the weights stay.
void stopInRound()
{
    parts::Part part = flatPart({}, 0, {}, true);
    EXPECT(part.system.has_value(), true);
    if (!part.system) {
        return;
    }
    param::System& system = *part.system;
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    const param::SystemContext context(part, false);
    NoValues regions;
    NoValues borders;
    int reports = 0;
    const std::function<void(float)> report = [&reports](float) { ++reports; };
    param::SolveStrategy strategy(context, system, regions, borders, progress, report);
    solver::ConstrainedDriver driver(system.normal.matrix(), system.normal.rightHandSide(), strategy, progress);
    driver.setIntegerUnknowns(system.integerUnknowns, system.integerPeriods);
    strategy.attach(driver);
    const std::vector<double> x(static_cast<std::size_t>(system.unknownCount()), 0.0);
    const std::vector<double> weights = system.weights;
    solver::DriverRequest request;
    request.x = &x;
    cancel = true;
    bool cancelled = false;
    try {
        strategy.onRequest(request);
    } catch (const api::StatusError& error) {
        cancelled = error.status() == ARM_CANCELLED;
    }
    EXPECT(cancelled, true);
    EXPECT(strategy.rounds(), 0);
    EXPECT(reports, 0);
    EXPECT(regions.calls + borders.calls, 0);
    EXPECT(request.systemChanged, false);
    EXPECT(system.weights == weights, true);
}

// Rounding of what is left: every free position of the integer list gets the
// whole number nearest to its value, half away from zero, in the order of the
// positions. The number goes to double as it is: one above 2^24 stays.
void finalRoundingValues()
{
    parts::Part part = flatPart({}, 0, {}, false);
    EXPECT(part.system.has_value(), true);
    if (!part.system) {
        return;
    }
    const param::SystemContext context(part, false);
    QuietProgress q;
    solver::CompactMatrix matrix;
    std::vector<double> rhs;
    NoAnswer callback;
    solver::ConstrainedDriver driver(matrix, rhs, callback, q.progress);
    driver.setIntegerUnknowns({2, 0, 1});
    const std::vector<double> x{-2.5, 0.4, 16777217.0};
    solver::DriverRequest request;
    request.x = &x;
    param::FinalRounding rounding;
    EXPECT(rounding.select(param::SolveState{context, *part.system, driver, q.progress}, request) ==
               param::RoundingStep::Outcome::Delivered,
           true);
    std::vector<std::pair<int, double>> selected;
    for (const solver::DriverRequest::Selection& chosen : request.selection) {
        selected.emplace_back(chosen.position, chosen.value);
    }
    EXPECT(selected, (std::vector<std::pair<int, double>>{{0, 16777217.0}, {1, -3.0}, {2, 0.0}}));
}

// The values of the variables after an elimination that failed, on the
// system of systemPenaltyRows: three conditions that were not reached stand
// empty in the list, the fourth is x44 = -2 x24, and the variables 0 and 1
// are the anchor. Every unknown has its value from the solution, x44 follows
// from x24, the anchor is 0; the empty conditions write nothing, so no value
// lands at another variable.
void valuesAfterFailedElimination()
{
    QuietProgress q;
    parts::Part part = flatPart({}, 0, {}, false);
    param::SystemSetup setup(part, false, q.progress);
    setup.constrainChains();
    std::vector<param::TermList> rows{
        {{40, 2.0}, {20, 1.0}}, {{30, 0x1p-34}, {50, 2.0}}, {{44, 2.0}, {24, 4.0}}, {{44, 4.0}, {10, 1.0}}};
    EXPECT(setup.eliminate(rows, param::EliminationMode::IsoLines, {}, false), false);
    EXPECT(setup.assemble(rows) == param::SystemSetup::Outcome::Built, true);
    const param::System& system = setup.system();
    EXPECT(system.variableCount, 60);
    EXPECT(system.unknownCount(), 57);
    int withoutPivot = 0;
    for (const param::Constraint& condition : system.constraints.list()) {
        withoutPivot += condition.pivot < 0 ? 1 : 0;
    }
    EXPECT(withoutPivot, 3);

    std::vector<double> x;
    std::vector<float> expected(60, 0.0f);
    for (int unknown = 0; unknown < system.unknownCount(); ++unknown) {
        x.push_back(100.0 + unknown);
        expected[static_cast<std::size_t>(system.unknownToVariable[static_cast<std::size_t>(unknown)])] =
            100.0f + static_cast<float>(unknown);
    }
    EXPECT(system.variableToUnknown[44], -1);
    EXPECT(expected[24], 122.0f);
    expected[44] = -244.0f;
    EXPECT(param::UvOutput::variables(system, x), expected);
}

// An elimination at the chain ends that fails, with the rows of
// systemPenaltyRows in place of the end rows on the open disc. Its rows and
// conditions speak of variables again, as after the iso lines; the set-up
// does not hand the rows that are left on (no penalty rows after the chain
// ends), so the system is assembled without them. The part counts as
// solved: the solution is finite and keeps the reached condition
// x44 = -2 x24, the rows that were left need not hold.
void chainEndsFailed()
{
    using param::Constraint;
    using param::TermList;
    using Rows = std::vector<TermList>;
    const double tiny = 0x1p-34;
    QuietProgress q;
    parts::Part part = flatPart({}, 0, {}, false);
    param::SystemSetup setup(part, false, q.progress);
    setup.constrainChains();
    Rows rows{{{40, 2.0}, {20, 1.0}}, {{30, tiny}, {50, 2.0}}, {{44, 2.0}, {24, 4.0}}, {{44, 4.0}, {10, 1.0}}};
    EXPECT(setup.eliminate(rows, param::EliminationMode::ChainEnds, {}, true), false);
    EXPECT(flat(rows),
           flat(Rows{{{20, 1.0}, {40, 2.0}}, {{10, 1.0}, {44, 4.0}}, {{30, tiny}, {50, 2.0}}, {{44, 2.0}, {24, 4.0}}}));
    EXPECT(flat(setup.system().constraints.list()),
           flat(std::vector<Constraint>{{0, {}}, {1, {}}, {-1, {}}, {-1, {}}, {-1, {}}, {44, {{24, -2.0}}}}));
    EXPECT(setup.assemble({}) == param::SystemSetup::Outcome::Built, true);
    EXPECT(setup.system().penaltyRowCount, 0);
    EXPECT(setup.system().rows.listCount(), 160);
    EXPECT(setup.system().unknownCount(), 57);
    part.system = setup.take();
    param::SolutionBuilder::compute(part, param::SolutionBuilder::Options{}, q.progress);
    EXPECT(part.solution.has_value(), true);
    if (!part.solution) {
        return;
    }
    const param::UvOutput::Coordinates& uv = part.solution->uv;
    bool finite = true;
    for (const std::array<float, 2>& p : uv) {
        finite = finite && std::isfinite(p[0]) && std::isfinite(p[1]);
    }
    EXPECT(finite, true);
    EXPECT(uv[22][0], -2.0f * uv[12][0]);
    EXPECT(uv[10][0] + 2.0f * uv[20][0] != 0.0f, true);
}

// Coordinates of a mirrored half. Two squares side by side, mirror plane
// x = 1: the triangles (0, 1, 4) and (0, 4, 3) on the left are computed,
// their mirror triangles (2, 4, 1) and (2, 5, 4) take their coordinates over
// with u negated, corner by corner. The triangle (5, 6, 4) lies in the
// mirrored half too, but point 6 has no mirror point and the triangle no
// mirror triangle: its corners keep what they have. Coordinates before: the
// grid u = 2 (x - 1), v = 4y on the left, 0 at the points 2 and 5 that only
// the mirrored half touches, (7, 9) at point 6.
void mirroredCoordinates()
{
    parts::Part part;
    std::vector<Vec3f> p{{0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}, {2.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f},
                         {1.0f, 1.0f, 0.0f}, {2.0f, 1.0f, 0.0f}, {2.0f, 2.0f, 0.0f}};
    std::vector<int> corners{0, 1, 4, 0, 4, 3, 2, 4, 1, 2, 5, 4, 5, 6, 4};
    part.mesh.swapPoints(p);
    part.mesh.editTopology().swapCorners(corners, FaceMode::Triangles);
    part.mesh.topologyChanged();
    field::CutGraph cuts;
    cuts.cutEdgeBits = BitArray(static_cast<std::size_t>(part.mesh.adjacency().edgeCount()));
    cuts.copy = field::CutCopy::cut(part.mesh, cuts.cutEdges, cuts.cutEdgeBits);
    cuts.mirrorHalf = BitArray(5);
    for (const std::size_t t : {2u, 3u, 4u}) {
        cuts.mirrorHalf.set(t);
    }
    cuts.mirrorHalfCount = 3;
    part.cutGraph = std::move(cuts);
    const param::SystemContext context(part, false);
    prep::SymmetryCache pairs = prep::SymmetryCache::fromPoints({2, 1, 0, 5, 4, 3, -1}, BitArray(7), std::nullopt);
    EXPECT(pairs.faceMirror(part.mesh), (std::vector<int>{2, 3, 0, 1, -1}));

    using Coordinates = param::UvOutput::Coordinates;
    Coordinates uv{{-2.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {-2.0f, 4.0f}, {0.0f, 4.0f}, {0.0f, 0.0f}, {7.0f, 9.0f}};
    param::UvOutput::mirror(context, pairs, uv);
    EXPECT(uv, (Coordinates{{-2.0f, 0.0f}, {0.0f, 0.0f}, {2.0f, 0.0f}, {-2.0f, 4.0f}, {0.0f, 4.0f}, {2.0f, 4.0f},
                            {7.0f, 9.0f}}));
}

// Regions and bridges on the flat grid with its border as crease curves, for
// the perfect grid u = 2x, v = 4y (a step along x is 2 long, one along y 4).
//  - first coordinate: the left and the right side are the sources. A point
//    in column x lies 2x from the left and 10 - 2x from the right side, so
//    columns 0 to 2 belong to the left and 3 to 5 to the right side, each
//    point reached along its row. The cheapest edges between the regions
//    are the five from column 2 to column 3 (4 + 4 + 2); the first of them,
//    2-3, is the bridge. Its ends lead back to points 0 and 5: variables 0
//    and 10 in the order of the ends, 10 grid lines apart,
//  - second coordinate: the lower and the upper side are the sources, each
//    point is reached along its column, and the bridge joins the two ends
//    of column 0 over edge 6-12: v of the points 0 and 24 (variables 1 and
//    49), 16 grid lines apart.
void regionsAndBridges()
{
    parts::Part part = flatPart({}, 0, {}, true);
    const param::SystemContext context(part, false);
    const Adjacency& copy = context.copy().adjacency();
    const param::UvOutput::Coordinates uv = param::UvOutput::coordinates(*part.system, shearedGrid(part, 2.0, 0.0));
    const param::RegionGrowth first(context, uv, 0);
    const param::RegionGrowth second(context, uv, 1);
    EXPECT(first.chainRegions(), 4);
    EXPECT(first.region(0) == first.region(24) && first.region(0) == first.region(14), true);
    EXPECT(first.region(5) == first.region(29) && first.region(5) == first.region(15), true);
    EXPECT(first.region(0) != first.region(5) && first.region(0) < 4 && first.region(5) < 4, true);
    EXPECT(first.distance(12), 0.0f);
    EXPECT(first.distance(14), 4.0f);
    EXPECT(first.distance(15), 4.0f);
    EXPECT(first.root(14), 12);
    EXPECT(first.root(15), 17);
    EXPECT(first.pathToRoot(14), (std::vector<int>{14, 13, 12}));
    EXPECT(first.onChain(copy.edgeBetween(0, 6)), true);
    EXPECT(first.onChain(copy.edgeBetween(0, 1)), false);
    EXPECT(second.onChain(copy.edgeBetween(0, 1)), true);
    EXPECT(second.onChain(copy.edgeBetween(0, 6)), false);
    EXPECT(second.region(1) == second.region(7) && second.region(25) == second.region(19), true);
    EXPECT(second.region(1) != second.region(25), true);
    EXPECT(second.distance(7), 4.0f);
    EXPECT(second.root(7), 1);
    EXPECT(second.root(19), 25);

    param::RegionBridges bridges;
    bridges.add(first, uv, copy);
    bridges.add(second, uv, copy);
    EXPECT(bridges.count(), 2);
    if (bridges.count() != 2) {
        return;
    }
    const param::Bridge& across = bridges.list()[0];
    EXPECT(across.coordinate == 0 && across.lowIsChain && across.highIsChain, true);
    EXPECT(across.edge, copy.edgeBetween(2, 3));
    EXPECT(copy.edgeVertices(across.edge) == (std::array<int, 2>{2, 3}), true);
    // The first variable belongs to the source of the first end of the edge.
    EXPECT(across.first, 0);
    EXPECT(across.second, 10);
    EXPECT(across.span, 10.0f);
    const param::Bridge& upwards = bridges.list()[1];
    EXPECT(upwards.coordinate == 1 && upwards.lowIsChain && upwards.highIsChain, true);
    EXPECT(upwards.edge, copy.edgeBetween(6, 12));
    EXPECT(upwards.first, 1);
    EXPECT(upwards.second, 49);
    EXPECT(upwards.span, 16.0f);
    EXPECT(bridges.bySpan(), (std::vector<int>{0, 1}));
}

// The additional system of a part as the rounding by regions sets it up for
// the coordinates `uv` of the copy points, with the rows of the bridges;
// `bridges` gets the bridges.
param::SoftSystem bridgeRowsAt(parts::Part& part, const param::UvOutput::Coordinates& uv,
                               param::RegionBridges& bridges)
{
    const param::SystemContext context(part, false);
    const param::System& system = *part.system;
    const Adjacency& copy = context.copy().adjacency();
    const std::array<param::RegionGrowth, 2> growth{param::RegionGrowth(context, uv, 0),
                                                    param::RegionGrowth(context, uv, 1)};
    for (const param::RegionGrowth& regions : growth) {
        bridges.add(regions, uv, copy);
    }
    param::SoftSystem rows(system.variableCount, bridges.count());
    param::BridgeRows::add(context, system.constraints, uv, growth, bridges, rows);
    return rows;
}

// The same for the unknowns `x`.
param::SoftSystem bridgeSystem(parts::Part& part, const std::vector<double>& x)
{
    param::RegionBridges bridges;
    return bridgeRowsAt(part, param::UvOutput::coordinates(*part.system, x), bridges);
}

// The perfect grid u = 2x, v = 4y at every copy point of a flat grid part.
param::UvOutput::Coordinates perfectGrid(const parts::Part& part)
{
    param::UvOutput::Coordinates uv;
    for (const int p : part.cutGraph->copy.pointOrigins()) {
        uv.push_back({2.0f * static_cast<float>(p % 6), 4.0f * static_cast<float>(p / 6)});
    }
    return uv;
}

// Index of the soft row with exactly these terms in this order, -1: none.
int softRowWith(const param::SoftSystem& rows, const param::TermList& terms)
{
    for (std::size_t i = 0; i < rows.softRows().size(); ++i) {
        if (flat(rows.softRows()[i]) == flat(terms)) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

// Rows of the bridges for the kinds of pairs the flat grid with border
// curves does not have, term by term, with the coordinates of the copy
// points set by hand on the perfect grid u = 2x, v = 4y. A row asks for
// "second - first = target": the term of the second variable (+1) comes
// first, then that of the first (-1), an anchored variable is left out;
// a core that exchanged the two sides of a bridge would turn the terms and
// the target round.
//  - two singular points: singular points 9, 14 and 15 without curves (9 is
//    the anchor); point 15 is moved to point 14 + (a, b). The bridge 14-15
//    of the first coordinate (u of 15 - u of 14, variables 30 and 28) is a
//    along and b across, that of the second (variables 31 and 29) b along
//    and a across. One soft row each with the present difference; the
//    weight is 0.01 when one fifth of the distance across is more than the
//    distance along, else a ramp from 0.01 by the ratio across / along
//    when that is more than 1, else 100 below 1, falling to 10 at 2.5, to
//    1 at 5 and 1 beyond. A distance across that is no number
//    leaves the weight by the distance along; one along that is no number
//    gives weight 1,
//  - a curve and a singular point: the border as curves, singular points 27
//    (the anchor; on the upper curve, a source of the first coordinate only)
//    and 14:
//     - point 14 at (0.5, 8): left side (u of point 12, variable 24) and
//       point 14 (variable 28) are 0.5 apart: a soft row "= 1" and a hard
//       row with the count of the bridge (column 62); point 27 at u = 6 and
//       the right side 4 apart: a soft row and a count row (column 65) in
//       which only the right side is left; 27 and 14 two singular points
//       mostly across; the other pairs only soft rows,
//     - point 14 at (3, 8): 3 from the left side, soft row and count row;
//       7 from the right side, a soft row only,
//     - point 14 at (5, 0.25): 0.25 above the lower side, but 1 beside the
//       nearest point of it: a soft row only,
//  - two curves that meet: the frames of the upper half of the grid turned
//    by a half turn, so that each side curve splits where the label of its
//    edges changes. Of the right side, the run 5-11-17 and the run 17-23-29
//    share point 17, which the later run (the lower one) takes; the edge
//    17-23 joins the two sources, the bridge has no edge and no row. Of
//    the left side, the edge 12-18 is left out of the curves: the two runs
//    are 0 apart along and 4 across, a soft row "= 0".
void bridgeRowClasses()
{
    using param::TermList;
    {
        parts::Part part = flatPart({}, 0, {9, 14, 15}, false);
        EXPECT(part.system.has_value(), true);
        if (!part.system) {
            return;
        }
        const float nan = std::numeric_limits<float>::quiet_NaN();
        struct Case {
            std::array<float, 2> offset;
            double weightAlongU;
            double weightAlongV;
        };
        const double low = static_cast<double>(0.01f);
        const std::array<Case, 6> cases{{
            {{0.5f, 0.0625f}, 100.0, low},
            {{2.0f, 0.25f}, 40.0, low},
            {{3.75f, 0.5f}, 5.5, low},
            {{6.0f, 1.0f}, 1.0, low},
            {{2.0f, 3.0f}, static_cast<double>(((3.0f / 2.0f - 0.2f) * 0.09f) / 0.8f + 0.01f),
             static_cast<double>(((3.0f - 2.5f) * -9.0f) / 2.5f + 10.0f)},
            {{2.0f, nan}, 40.0, 1.0},
        }};
        for (const Case& c : cases) {
            param::UvOutput::Coordinates uv = perfectGrid(part);
            uv[15] = {uv[14][0] + c.offset[0], uv[14][1] + c.offset[1]};
            param::RegionBridges bridges;
            const param::SoftSystem rows = bridgeRowsAt(part, uv, bridges);
            EXPECT(rows.hardRows().empty(), true);
            const int u = softRowWith(rows, TermList{{30, 1.0}, {28, -1.0}});
            const int v = softRowWith(rows, TermList{{31, 1.0}, {29, -1.0}});
            EXPECT(u >= 0 && v >= 0, true);
            if (u < 0 || v < 0) {
                continue;
            }
            EXPECT(rows.softTargets()[static_cast<std::size_t>(u)], static_cast<double>(c.offset[0]));
            EXPECT(rows.softWeights()[static_cast<std::size_t>(u)], c.weightAlongU);
            const double targetV = rows.softTargets()[static_cast<std::size_t>(v)];
            EXPECT(std::isnan(c.offset[1]) ? std::isnan(targetV) : targetV == static_cast<double>(c.offset[1]), true);
            EXPECT(rows.softWeights()[static_cast<std::size_t>(v)], c.weightAlongV);
        }
    }
    {
        parts::Part part = flatPart({}, 0, {27, 14}, true);
        EXPECT(part.system.has_value(), true);
        if (!part.system) {
            return;
        }
        auto rowsWith = [&part](std::array<float, 2> place) {
            param::UvOutput::Coordinates uv = perfectGrid(part);
            uv[14] = place;
            param::RegionBridges bridges;
            return bridgeRowsAt(part, uv, bridges);
        };
        {
            const param::SoftSystem rows = rowsWith({0.5f, 8.0f});
            EXPECT(flat(rows.softRows()), flat(std::vector<TermList>{{{28, 1.0}, {24, -1.0}},
                                                                     {{58, 1.0}},
                                                                     {{28, -1.0}},
                                                                     {{48, -1.0}},
                                                                     {{29, 1.0}, {3, -1.0}},
                                                                     {{53, 1.0}, {29, -1.0}},
                                                                     {{34, 1.0}, {28, -1.0}},
                                                                     {{10, 1.0}, {0, -1.0}},
                                                                     {{59, 1.0}, {11, -1.0}}}));
            EXPECT(rows.softTargets(), (std::vector<double>{1.0, 4.0, 5.5, 6.0, 8.0, 8.0, 9.5, 10.0, 16.0}));
            const double askew = static_cast<double>(((8.0f / 5.5f - 0.2f) * 0.09f) / 0.8f + 0.01f);
            EXPECT(rows.softWeights(), (std::vector<double>{1.0, 1.0, askew, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}));
            EXPECT(flat(rows.hardRows()),
                   flat(std::vector<TermList>{{{24, -1.0}, {28, 1.0}, {62, -1.0}}, {{58, 1.0}, {65, -1.0}}}));
            EXPECT(rows.hardTargets(), (std::vector<double>{0.0, 0.0}));
        }
        {
            const param::SoftSystem rows = rowsWith({3.0f, 8.0f});
            const int left = softRowWith(rows, TermList{{28, 1.0}, {24, -1.0}});
            const int right = softRowWith(rows, TermList{{34, 1.0}, {28, -1.0}});
            EXPECT(left >= 0 && right >= 0, true);
            if (left >= 0 && right >= 0) {
                EXPECT(rows.softTargets()[static_cast<std::size_t>(left)], 3.0);
                EXPECT(rows.softTargets()[static_cast<std::size_t>(right)], 7.0);
            }
            EXPECT(rows.hardRows().empty() ? std::vector<std::uint64_t>{} : flat(rows.hardRows().front()),
                   flat(TermList{{24, -1.0}, {28, 1.0}, {62, -1.0}}));
            EXPECT(rows.hardRows().size(), std::size_t{2});
        }
        {
            const param::SoftSystem rows = rowsWith({5.0f, 0.25f});
            const int above = softRowWith(rows, TermList{{29, 1.0}, {5, -1.0}});
            EXPECT(above >= 0, true);
            if (above >= 0) {
                EXPECT(rows.softTargets()[static_cast<std::size_t>(above)], 0.25);
                EXPECT(rows.softWeights()[static_cast<std::size_t>(above)], 1.0);
            }
            EXPECT(flat(rows.hardRows()), flat(std::vector<TermList>{{{58, 1.0}, {65, -1.0}}}));
        }
    }
    {
        parts::Part part = flatPart({}, 0, {}, true);
        for (std::size_t t = 20; t < 40; ++t) {
            part.crossField->frameE1[t] = part.crossField->frameE1[t] * -1.0f;
            part.crossField->frameE2[t] = part.crossField->frameE2[t] * -1.0f;
        }
        field::CurveAlignment::orient(part.mesh, *part.crossField, *part.cutGraph);
        QuietProgress q;
        param::ChainBuilder::compute(part, param::ChainBuilder::Options{}, q.progress);
        param::SystemBuilder::compute(part, param::SystemBuilder::Options{}, q.progress);
        EXPECT(part.system.has_value(), true);
        if (!part.system) {
            return;
        }
        const param::SystemContext context(part, false);
        const Adjacency& copy = context.copy().adjacency();
        const Csr& curves = part.chains->curves.edges;
        EXPECT(part.chains->curves.count(), 6);
        if (part.chains->curves.count() != 6) {
            return;
        }
        EXPECT(v(curves.row(0)), edgesBetween(copy, {{23, 29}, {17, 23}}));
        EXPECT(v(curves.row(1)), edgesBetween(copy, {{11, 17}, {5, 11}}));
        EXPECT(v(curves.row(2)), edgesBetween(copy, {{18, 24}}));
        EXPECT(v(curves.row(3)), edgesBetween(copy, {{6, 12}, {0, 6}}));
        const param::UvOutput::Coordinates uv = perfectGrid(part);
        const param::RegionGrowth first(context, uv, 0);
        EXPECT(first.region(17), 1);
        EXPECT(first.region(23), 0);
        param::RegionBridges bridges;
        const param::SoftSystem rows = bridgeRowsAt(part, uv, bridges);
        int meeting = 0;
        for (const param::Bridge& bridge : bridges.list()) {
            if (bridge.coordinate == 0 && bridge.low == 0 && bridge.high == 1) {
                ++meeting;
                EXPECT(bridge.edge == -1 && bridge.first == -1 && bridge.second == -1, true);
            }
        }
        EXPECT(meeting, 1);
        EXPECT(flat(rows.softRows()), flat(std::vector<TermList>{{{36, 1.0}, {24, -1.0}},
                                                                 {{10, 1.0}},
                                                                 {{46, 1.0}, {24, -1.0}},
                                                                 {{46, 1.0}, {36, -1.0}},
                                                                 {{49, 1.0}}}));
        EXPECT(rows.softTargets(), (std::vector<double>{0.0, 10.0, 10.0, 10.0, 16.0}));
        EXPECT(rows.hardRows().empty(), true);
    }
}

// The field sum that gives a bridge its side divides by the size along the
// coordinate of the bridge. Flat grid with border curves, v = 0.75 y (the
// upper side 3 above the lower one: a soft row and a count row); in the
// upper half (rows 2 and 3) the second combed direction is turned round and
// its size is 0.125, below it stays 0.25. The bridge of the second
// coordinate runs up column 0 from point 0 to point 24: the field counts
// 4 + 4 steps in the lower half and -8 - 8 in the upper, -8 in all, against
// the difference of +3 - so the difference counts as -3. Divided by the size
// of the first direction (0.5 everywhere) the sum would be 0 and the
// difference would stay +3.
void bridgeFieldSizes()
{
    using param::TermList;
    parts::Part part = flatPart({}, 0, {}, true);
    for (std::size_t t = 20; t < 40; ++t) {
        part.cutGraph->combedV[t] = {-part.cutGraph->combedV[t][0], -part.cutGraph->combedV[t][1]};
        part.sizeField->sizeV[t] = 0.125f;
    }
    param::UvOutput::Coordinates uv = perfectGrid(part);
    for (std::array<float, 2>& point : uv) {
        point[1] = point[1] * 0.1875f;
    }
    param::RegionBridges bridges;
    const param::SoftSystem rows = bridgeRowsAt(part, uv, bridges);
    EXPECT(bridges.count(), 2);
    if (bridges.count() != 2) {
        return;
    }
    const param::Bridge& upwards = bridges.list()[1];
    EXPECT(upwards.edge, part.cutGraph->copy.adjacency().edgeBetween(6, 12));
    EXPECT(upwards.first == 1 && upwards.second == 49, true);
    EXPECT(flat(rows.softRows()), flat(std::vector<TermList>{{{49, 1.0}}, {{10, 1.0}}}));
    EXPECT(rows.softTargets(), (std::vector<double>{-3.0, 10.0}));
    EXPECT(flat(rows.hardRows()), flat(std::vector<TermList>{{{49, -1.0}, {62, -1.0}}}));
}

// Rows of the additional system on the flat grid with border curves (60
// variables: constant column 60, bridge columns 61 and 62; variable 0 is
// anchored, variable 10 is u at point 5 on the right side). The bridge of
// the first coordinate runs along the lower border, where the field counts
// 2 grid steps per edge in the direction of x.
//  - perfect grid (right side at 10, upper side at 16): both bridges join
//    two curves more than 5 apart: one soft row each with weight 1. Nothing
//    to eliminate; the two references of the sides are the unknowns, each
//    tied to its present value with the small weight,
//  - right side at 3: between 1.35 and 5 and straight across: the soft row
//    and a hard row "u(5) = count of bridge 0". The elimination writes the
//    reference of the right side as that count, which becomes the unknown,
//    tied to 1,
//  - right side at 1: closer than 1.35: only a hard row "u(5) = 1". Its
//    target goes to the constant column: the reference of the right side is
//    minus the constant column, which stands for -1.
void softSystemRows()
{
    using param::TermList;
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    const double tie = static_cast<double>(1e-4f);
    {
        parts::Part part = flatPart({}, 0, {}, true);
        const std::vector<double> x = shearedGrid(part, 2.0, 0.0);
        param::SoftSystem rows = bridgeSystem(part, x);
        EXPECT(rows.constantColumn(), 60);
        EXPECT(rows.bridgeColumn(1), 62);
        EXPECT(rows.columnCount(), 63);
        EXPECT(rows.hardRows().empty(), true);
        EXPECT(rows.softWeights(), (std::vector<double>{1.0, 1.0}));
        const bool two = rows.softRows().size() == 2 && rows.softRows()[0].size() == 1;
        EXPECT(two, true);
        if (!two) {
            return;
        }
        // Term of the second variable +1 (the anchored first one is left
        // out), target its present difference to the first.
        EXPECT(flat(rows.softRows()[0]), flat(TermList{{10, 1.0}}));
        EXPECT(flat(rows.softRows()[1]), flat(TermList{{49, 1.0}}));
        EXPECT(rows.softTargets(), (std::vector<double>{10.0, 16.0}));

        const param::ReducedSoftSystem reduced = std::move(rows).reduce(*part.system, {}, x, progress);
        EXPECT(reduced.conditions.count(), 0);
        EXPECT(reduced.unknownCount(), 2);
        EXPECT(reduced.rows.listCount(), 4);
        EXPECT(reduced.rows.otherDimension(), 2);
        EXPECT(reduced.weights, (std::vector<double>{1.0, 1.0, tie, tie}));
        if (reduced.unknownCount() != 2 || reduced.rows.listCount() != 4) {
            return;
        }
        EXPECT(reduced.columnOf[0] % 2 == 0 && reduced.columnOf[1] % 2 == 1, true);
        EXPECT(reduced.unknownOf[static_cast<std::size_t>(reduced.columnOf[1])], 1);
        const std::array<double, 2> span{10.0, 16.0};
        for (int row = 0; row < 2; ++row) {
            const TermList& terms = reduced.rows.list(row);
            const bool single = terms.size() == 1 && terms[0].index == row;
            EXPECT(single, true);
            if (single) {
                EXPECT(terms[0].value * reduced.targets[static_cast<std::size_t>(row)],
                       span[static_cast<std::size_t>(row)]);
            }
        }
        for (int row = 2; row < 4; ++row) {
            const TermList& terms = reduced.rows.list(row);
            const bool single = terms.size() == 1 && terms[0].value == 1.0;
            EXPECT(single, true);
            if (single) {
                EXPECT(reduced.targets[static_cast<std::size_t>(row)],
                       span[static_cast<std::size_t>(terms[0].index)]);
            }
        }
    }
    {
        parts::Part part = flatPart({}, 0, {}, true);
        const std::vector<double> x = shearedGrid(part, 0.6, 0.0);
        param::SoftSystem rows = bridgeSystem(part, x);
        const bool two = rows.softRows().size() == 2 && rows.softRows()[0].size() == 1;
        EXPECT(two, true);
        if (!two) {
            return;
        }
        EXPECT(flat(rows.softRows()[0]), flat(TermList{{10, 1.0}}));
        EXPECT(rows.softTargets()[0], 3.0);
        EXPECT(flat(rows.hardRows()), flat(std::vector<TermList>{{{10, 1.0}, {61, -1.0}}}));
        EXPECT(rows.hardTargets(), (std::vector<double>{0.0}));

        const param::ReducedSoftSystem reduced = std::move(rows).reduce(*part.system, {}, x, progress);
        EXPECT(reduced.conditions.count(), 1);
        EXPECT(reduced.unknownCount(), 2);
        EXPECT(reduced.rows.listCount(), 4);
        if (reduced.conditions.count() != 1 || reduced.unknownCount() != 2 || reduced.rows.listCount() != 4) {
            return;
        }
        const param::Constraint& condition = reduced.conditions.list()[0];
        EXPECT(flat(condition.terms), flat(TermList{{61, 1.0}}));
        EXPECT(condition.pivot % 2 == 0 && part.system->variableToUnknown[static_cast<std::size_t>(condition.pivot)] >= 0,
               true);
        EXPECT(reduced.conditions.ofVariable(condition.pivot), 0);
        EXPECT(reduced.unknownOf[static_cast<std::size_t>(condition.pivot)], -1);
        EXPECT(reduced.columnOf[0], 61);
        const bool counted = reduced.rows.list(0).size() == 1 && reduced.rows.list(0)[0].index == 0;
        EXPECT(counted, true);
        if (counted) {
            EXPECT(reduced.rows.list(0)[0].value * reduced.targets[0], 3.0);
        }
        // The rows that tie: the reference of the upper side to 16, the count to 1.
        EXPECT(flat(reduced.rows.list(2)), flat(TermList{{1, 1.0}}));
        EXPECT(flat(reduced.rows.list(3)), flat(TermList{{0, 1.0}}));
        EXPECT((std::vector<double>(reduced.targets.begin() + 2, reduced.targets.end())),
               (std::vector<double>{16.0, 1.0}));
    }
    {
        parts::Part part = flatPart({}, 0, {}, true);
        const std::vector<double> x = shearedGrid(part, 0.2, 0.0);
        param::SoftSystem rows = bridgeSystem(part, x);
        EXPECT(rows.softRows().size(), std::size_t{1});
        const bool single = rows.hardRows().size() == 1 && rows.hardRows()[0].size() == 1;
        EXPECT(single, true);
        if (!single) {
            return;
        }
        EXPECT(flat(rows.hardRows()[0]), flat(TermList{{10, 1.0}}));
        EXPECT(rows.hardTargets(), (std::vector<double>{1.0}));

        const param::ReducedSoftSystem reduced = std::move(rows).reduce(*part.system, {}, x, progress);
        EXPECT(reduced.conditions.count(), 1);
        EXPECT(reduced.unknownCount(), 1);
        if (reduced.conditions.count() == 1) {
            EXPECT(flat(reduced.conditions.list()[0].terms), flat(TermList{{60, -1.0}}));
            EXPECT(reduced.conditions.list()[0].pivot % 2, 0);
        }
    }
}

// What the rounding by regions selects for a part and the unknowns `x`,
// asked outside a solve (no position is fixed).
struct RegionAnswer {
    param::RoundingStep::Outcome outcome = param::RoundingStep::Outcome::Failed;
    std::vector<solver::DriverRequest::Selection> selection;
};

RegionAnswer askRegions(parts::Part& part, const std::vector<double>& x)
{
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    param::System& system = *part.system;
    const param::SystemContext context(part, false);
    NoAnswer idle;
    solver::ConstrainedDriver driver(system.normal.matrix(), system.normal.rightHandSide(), idle, progress);
    driver.setIntegerUnknowns(system.integerUnknowns, system.integerPeriods);
    solver::DriverRequest request;
    request.x = &x;
    param::RegionRounding regions;
    RegionAnswer answer;
    answer.outcome = regions.select(param::SolveState{context, system, driver, progress}, request);
    answer.selection = request.selection;
    return answer;
}

// A solve that lets the driver round one unknown itself and asks the
// regions in the round after it.
class RegionsSecond final : public solver::DriverCallback {
public:
    RegionsSecond(const param::SystemContext& context, const param::System& system, const run::Progress& progress)
        : context_(context), system_(system), progress_(progress)
    {
    }
    void attach(const solver::ConstrainedDriver& driver) { driver_ = &driver; }

    void onRequest(solver::DriverRequest& request) override
    {
        if (request.round == 0) {
            return;
        }
        fixedBefore = driver_->fixedPositions().testBounded(0) ? 0 : 1;
        param::RegionRounding regions;
        outcome = regions.select(param::SolveState{context_, system_, *driver_, progress_}, request);
        selection = request.selection;
        request.answer = solver::DriverRequest::Answer::Delivered;
    }

    param::RoundingStep::Outcome outcome = param::RoundingStep::Outcome::Failed;
    std::vector<solver::DriverRequest::Selection> selection;
    int fixedBefore = -1;

private:
    const param::SystemContext& context_;
    const param::System& system_;
    const run::Progress& progress_;
    const solver::ConstrainedDriver* driver_ = nullptr;
};

// Values the rounding by regions hands back on the flat grid with border
// curves; the two unknowns to round are u of the right side and v of the
// upper side (always 16 here).
//  - right side at 10: taken over from the additional system, where the
//    batch rounding makes it whole,
//  - right side at 3: the count of the bridge is rounded to 3 and goes back
//    through the condition of the hard row,
//  - right side at 1: the hard row gives 1 through the constant column,
//  - right side at -3: the grid runs against the field; the difference
//    counts as positive in the direction of the field, the count is 3 and
//    the right side goes to +3.
// A position the solve has fixed before is left out.
void regionSelection()
{
    using Outcome = param::RoundingStep::Outcome;
    const std::array<std::array<double, 2>, 4> cases{{{2.0, 10.0}, {0.6, 3.0}, {0.2, 1.0}, {-0.6, 3.0}}};
    for (const std::array<double, 2>& scaleAndValue : cases) {
        parts::Part part = flatPart({}, 0, {}, true);
        const RegionAnswer answer = askRegions(part, shearedGrid(part, scaleAndValue[0], 0.0));
        EXPECT(answer.outcome == Outcome::Delivered, true);
        EXPECT(answer.selection.size(), std::size_t{2});
        if (answer.selection.size() != 2) {
            continue;
        }
        EXPECT(answer.selection[0].position + answer.selection[1].position, 1);
        for (const solver::DriverRequest::Selection& chosen : answer.selection) {
            const int unknown = part.system->integerUnknowns[static_cast<std::size_t>(chosen.position)];
            const int variable = part.system->unknownToVariable[static_cast<std::size_t>(unknown)];
            EXPECT(chosen.value, variable % 2 == 0 ? scaleAndValue[1] : 16.0);
        }
    }
    // With a grid period of 2 for u of the right side: right side at 11.4,
    // an unknown of the additional system, is rounded by the inner driver to
    // the period - 12; right side at 3 follows from the count of its bridge
    // (3) and is rounded to the period on its way back - 4. v of the upper
    // side keeps period 1.
    for (const std::array<double, 2>& scaleAndValue : {std::array<double, 2>{2.28, 12.0}, {0.6, 4.0}}) {
        parts::Part part = flatPart({}, 0, {}, true);
        param::System& system = *part.system;
        system.periods[46] = 2;
        system.integerPeriods = {2.0, 1.0};
        const RegionAnswer answer = askRegions(part, shearedGrid(part, scaleAndValue[0], 0.0));
        EXPECT(answer.outcome == Outcome::Delivered, true);
        EXPECT(answer.selection.size(), std::size_t{2});
        for (const solver::DriverRequest::Selection& chosen : answer.selection) {
            const int unknown = system.integerUnknowns[static_cast<std::size_t>(chosen.position)];
            const int variable = system.unknownToVariable[static_cast<std::size_t>(unknown)];
            EXPECT(chosen.value, variable == 46 ? scaleAndValue[1] : 16.0);
        }
    }
    {
        parts::Part part = flatPart({}, 0, {}, true);
        std::atomic<bool> cancel{false};
        run::Progress progress(nullptr, nullptr, cancel);
        param::System& system = *part.system;
        const param::SystemContext context(part, false);
        RegionsSecond rounds(context, system, progress);
        solver::ConstrainedDriver driver(system.normal.matrix(), system.normal.rightHandSide(), rounds, progress);
        driver.setIntegerUnknowns(system.integerUnknowns, system.integerPeriods);
        rounds.attach(driver);
        std::vector<double> x;
        EXPECT(driver.run(x), true);
        EXPECT(rounds.outcome == Outcome::Delivered, true);
        EXPECT(rounds.selection.size(), std::size_t{1});
        if (rounds.selection.size() == 1) {
            EXPECT(rounds.selection[0].position, 1 - rounds.fixedBefore);
            const int unknown = system.integerUnknowns[static_cast<std::size_t>(rounds.selection[0].position)];
            EXPECT(rounds.selection[0].value, system.unknownToVariable[static_cast<std::size_t>(unknown)] % 2 == 0 ? 10.0 : 16.0);
        }
    }
}

// The rounding by regions of the core, with what it answered in every round.
// Its first `failing` calls fail without asking the core.
class RecordedRegions final : public param::RoundingStep {
public:
    explicit RecordedRegions(int failing = 0) : failing_(failing) {}

    Outcome select(const param::SolveState& state, solver::DriverRequest& request) override
    {
        const bool fail = static_cast<int>(outcomes.size()) < failing_;
        outcomes.push_back(fail ? Outcome::Failed : regions_.select(state, request));
        selected.push_back(request.selection.size());
        return outcomes.back();
    }
    std::vector<Outcome> outcomes;
    std::vector<std::size_t> selected;

private:
    param::RegionRounding regions_;
    int failing_;
};

// Order of the rounds with the rounding by regions of the core.
//  - flat grid with border curves: reweighting, then the regions deliver
//    both unknowns; nothing is left for a further round: two rounds,
//  - cut chain from the border to a singular point without curves: a single
//    region per coordinate, so no bridge and an additional system without
//    unknowns. The regions fail in every round and the driver rounds one
//    unknown each time: three rounds for two unknowns.
void regionRounds()
{
    using Outcome = param::RoundingStep::Outcome;
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    {
        parts::Part part = flatPart({}, 0, {}, true);
        RecordedRegions regions;
        FirstPosition borders;
        param::SolutionBuilder::compute(part, param::SolutionBuilder::Options{}, progress,
                                        param::SolutionBuilder::Steps{regions, borders});
        EXPECT(part.solution.has_value(), true);
        if (!part.solution) {
            return;
        }
        EXPECT(regions.outcomes == (std::vector<Outcome>{Outcome::Delivered}), true);
        EXPECT(regions.selected, (std::vector<std::size_t>{2}));
        EXPECT(borders.calls, 0);
        EXPECT(part.solution->rounds, 2);
        EXPECT(gridMisfit(part) < 1e-4, true);
        EXPECT(notWhole(part), 0);
    }
    {
        parts::Part part = flatPart({{12, 13}, {13, 14}, {14, 15}}, 0, {15}, false);
        RecordedRegions regions;
        FirstPosition borders;
        param::SolutionBuilder::compute(part, param::SolutionBuilder::Options{}, progress,
                                        param::SolutionBuilder::Steps{regions, borders});
        EXPECT(part.solution.has_value(), true);
        if (!part.solution) {
            return;
        }
        EXPECT(regions.outcomes == (std::vector<Outcome>{Outcome::Failed, Outcome::Failed}), true);
        EXPECT(regions.selected, (std::vector<std::size_t>{0, 0}));
        EXPECT(part.solution->rounds, 3);
        EXPECT(notWhole(part), 0);
    }
}

// The rounding by regions with nothing left to deliver, on the flat grid
// with border curves and a first direction of size 5/11 (u = 2.2 x, the
// right side at 11). The integer list gets u of point 1 (2.2, variable 2)
// in front of u of the right side and v of the upper side; only these two
// are known to the regions.
//  - the regions fail in their first two calls: the driver rounds the two
//    unknowns nearest to a whole number, u of the right side and v of the
//    upper side. In the third call the regions are solved, but both their
//    positions are fixed: nothing is selected, the answer says so, and the
//    final rounding takes u of point 1 in the same round - four rounds,
//    the corner opposite the anchor at (11, 16), point 1 at u = 2,
//  - with u of point 1 as the only position the regions are solved but
//    have a value for no position: the empty selection ends the solve
//    without a solution.
void regionsAfterFixedPositions()
{
    using Outcome = param::RoundingStep::Outcome;
    QuietProgress q;
    auto stretched = [&q](bool withRegions) {
        parts::Part part = flatPart({}, 0, {}, true);
        part.sizeField->sizeU.assign(part.sizeField->sizeU.size(), 5.0f / 11.0f);
        param::SystemBuilder::compute(part, param::SystemBuilder::Options{}, q.progress);
        if (part.system) {
            param::System& system = *part.system;
            const int pointOne = system.variableToUnknown[2];
            if (withRegions) {
                system.integerUnknowns.insert(system.integerUnknowns.begin(), pointOne);
            } else {
                system.integerUnknowns = {pointOne};
            }
        }
        return part;
    };
    {
        parts::Part part = stretched(true);
        EXPECT(part.system.has_value(), true);
        if (!part.system) {
            return;
        }
        std::vector<int> whole;
        for (const int unknown : part.system->integerUnknowns) {
            whole.push_back(part.system->unknownToVariable[static_cast<std::size_t>(unknown)]);
        }
        EXPECT(whole, (std::vector<int>{2, 46, 57}));
        RecordedRegions regions(2);
        FirstPosition borders;
        param::SolutionBuilder::compute(part, param::SolutionBuilder::Options{}, q.progress,
                                        param::SolutionBuilder::Steps{regions, borders});
        EXPECT(part.solution.has_value(), true);
        if (!part.solution) {
            return;
        }
        EXPECT(regions.outcomes == (std::vector<Outcome>{Outcome::Failed, Outcome::Failed, Outcome::AlreadyFixed}),
               true);
        EXPECT(regions.selected, (std::vector<std::size_t>{0, 0, 0}));
        EXPECT(borders.calls, 0);
        EXPECT(part.solution->rounds, 4);
        EXPECT(part.solution->uv[29] == (std::array<float, 2>{11.0f, 16.0f}), true);
        EXPECT(part.solution->uv[1][0], 2.0f);
    }
    {
        parts::Part part = stretched(false);
        RecordedRegions regions;
        FirstPosition borders;
        param::SolutionBuilder::compute(part, param::SolutionBuilder::Options{}, q.progress,
                                        param::SolutionBuilder::Steps{regions, borders});
        EXPECT(regions.outcomes == (std::vector<Outcome>{Outcome::Delivered}), true);
        EXPECT(regions.selected, (std::vector<std::size_t>{0}));
        EXPECT(part.solution.has_value(), false);
    }
}

// Hard rows across a cut chain, with sides and coordinates set by hand:
// side A runs over copy points 3, 4, 5 and side B over 9, 10, 11; point 3
// lies at (1, 2), point 9 is moved by `offset` and point 5 by `run`. The
// row ties the coordinate of the smaller part of the offset at points 3
// and 9 (variables 6 and 18, or 7 and 19).
//  - smaller part below 1: offset above 4 and run above 6, or offset above 6
//    and run above 4, or both above 2 and the smaller part below 0.02 of
//    0.75 offset + 0.25 run,
//  - smaller part from 1 to below 2: offset above 14, run above 10 and the
//    smaller part below 0.08 of that mix,
//  - no row for a chain with a turn (a whole turn is none) and for a chain
//    along a disc cut.
void seamRowRules()
{
    using param::TermList;
    parts::Part part = flatPart({{12, 13}, {13, 14}, {14, 15}}, 0, {15}, false);
    const param::SystemContext context(part, false);
    param::Chains& chains = *part.chains;
    EXPECT(chains.count(), 1);
    if (chains.count() != 1) {
        return;
    }
    chains.sideA[0] = {3, 4, 5};
    chains.sideB[0] = {9, 10, 11};
    auto rowsFor = [&](float turns, std::array<float, 2> offset, std::array<float, 2> run) {
        chains.rotation[0] = turns;
        param::UvOutput::Coordinates uv(static_cast<std::size_t>(part.system->copyPointCount),
                                        std::array<float, 2>{1.0f, 2.0f});
        uv[9] = {1.0f + offset[0], 2.0f + offset[1]};
        uv[5] = {1.0f + run[0], 2.0f + run[1]};
        param::SoftSystem rows(part.system->variableCount, 0);
        param::SeamRows::add(context, uv, rows);
        EXPECT(rows.softRows().empty(), true);
        EXPECT(rows.hardTargets(), std::vector<double>(rows.hardRows().size(), 0.0));
        return flat(rows.hardRows());
    };
    const std::vector<std::uint64_t> none = flat(std::vector<TermList>{});
    const std::vector<std::uint64_t> alongFirst = flat(std::vector<TermList>{{{6, 1.0}, {18, -1.0}}});
    const std::vector<std::uint64_t> alongSecond = flat(std::vector<TermList>{{{7, 1.0}, {19, -1.0}}});
    EXPECT(rowsFor(0.0f, {5.0f, 0.5f}, {0.0f, 7.0f}), alongSecond);
    EXPECT(rowsFor(0.0f, {5.0f, 0.5f}, {0.0f, 5.0f}), none);
    EXPECT(rowsFor(0.0f, {7.0f, 0.5f}, {0.0f, 5.0f}), alongSecond);
    EXPECT(rowsFor(0.0f, {3.0f, 0.03125f}, {0.0f, 3.0f}), alongSecond);
    EXPECT(rowsFor(0.0f, {3.0f, 0.5f}, {0.0f, 3.0f}), none);
    EXPECT(rowsFor(0.0f, {15.0f, 1.0f}, {0.0f, 11.0f}), alongSecond);
    EXPECT(rowsFor(0.0f, {15.0f, 1.5f}, {0.0f, 11.0f}), none);
    EXPECT(rowsFor(0.0f, {15.0f, 2.0f}, {0.0f, 40.0f}), none);
    EXPECT(rowsFor(0.0f, {0.5f, 5.0f}, {7.0f, 0.0f}), alongFirst);
    EXPECT(rowsFor(0.0f, {0.5f, 5.0f}, {0.0f, 7.0f}), none);
    EXPECT(rowsFor(1.0f, {5.0f, 0.5f}, {0.0f, 7.0f}), none);
    EXPECT(rowsFor(-1.0f, {5.0f, 0.5f}, {0.0f, 7.0f}), none);
    EXPECT(rowsFor(4.0f, {5.0f, 0.5f}, {0.0f, 7.0f}), alongSecond);

    BitArray diskCuts(static_cast<std::size_t>(part.mesh.adjacency().edgeCount()));
    diskCuts.set(static_cast<std::size_t>(chains.edges.row(0)[0]));
    part.cutGraph->diskCutBits = diskCuts;
    EXPECT(rowsFor(0.0f, {5.0f, 0.5f}, {0.0f, 7.0f}), none);
}

} // namespace

namespace {

// Freezes the borders of a part: the points in `kept` are the border points
// that stay, every other point has no mark.
void freezeBorders(parts::Part& part, std::initializer_list<int> kept)
{
    std::vector<FreezeMark> marks(static_cast<std::size_t>(part.mesh.pointCount()));
    for (const int p : kept) {
        marks[static_cast<std::size_t>(p)].value = 1;
    }
    part.mesh.maps().set(std::make_unique<TypedMap<FreezeMark>>(std::string(map_names::freezeMarks),
                                                                ElementKind::Point, std::move(marks)));
    part.bordersFrozen = true;
}

// Asks the rounding along frozen borders and keeps what it selected, as
// (variable, value) in the order of the selection.
class RecordedBorders final : public param::RoundingStep {
public:
    Outcome select(const param::SolveState& state, solver::DriverRequest& request) override
    {
        ++calls;
        const std::size_t before = request.selection.size();
        outcome = step.select(state, request);
        for (std::size_t k = before; k < request.selection.size(); ++k) {
            const auto position = static_cast<std::size_t>(request.selection[k].position);
            const auto unknown = static_cast<std::size_t>(state.driver.integerUnknowns()[position]);
            picks.emplace_back(state.system.unknownToVariable[unknown], request.selection[k].value);
            alreadyFixed += state.driver.fixedPositions().testBounded(position) ? 1 : 0;
        }
        return outcome;
    }

    param::BorderRounding step;
    Outcome outcome = Outcome::Failed;
    int calls = 0;
    int alreadyFixed = 0;
    std::vector<std::pair<int, double>> picks;
};

// What the rounding along frozen borders selects, with system, marks and
// solution set by hand. The mesh is the flat grid without the quad between
// the points 8, 9, 15 and 14 and with a cut from border point 12 to point
// 13: border point 12 has two copies (12 below the cut, 30 above it) that no
// copy edge joins; copy 30 is joined to point 18, copy 12 to point 6. Two
// border loops: the outer one (0 1 2 3 4 5 11 17 23 29 28 27 26 25 24 18 12
// 6) and the one around the hole (9 8 14 15). Kept are the points 0 to 5,
// 18, 12 and 6 of the outer loop and point 14 of the hole. Variable 2w is u
// and 2w + 1 is v of copy point w, unknown = variable; every coordinate is a
// coordinate of a frozen border unless noted.
//  u along the outer loop, in list order:
//   - point 0: eliminated without terms, whole number 0,
//   - points 1, 2, 3 at 2.25, 2.375, 1.625 all round to 2: of 1 and 2 the
//     one farther away (2) keeps its value; 2 and 3 are equally far, so
//     again 2 - point 3 is rounded,
//   - point 4 at 4.5 rounds to 5, but its unknown is not in the integer
//     list: left out,
//   - point 5 at -0.5 rounds to -1,
//   - point 18 (a node, not a frozen coordinate) at 7 and copy 12 at 7.25
//     land on 7 but are not joined: both rounded,
//   - copy 30: eliminated, 1 * u(12) - 0.125 * u(1) = 6.96875, lands on 7
//     as well, not joined to copy 12 or point 6, and has no unknown,
//   - point 6 at 0.25 and point 0 (last and first of the list) land on 0:
//     point 0 has no unknown, so point 6 keeps its value.
//  v along the outer loop:
//   - point 0: eliminated, 0; point 1 at 0.125: point 0 has no unknown, so
//     point 1 keeps its value,
//   - point 2 at -0.25 (rounds to 0) is farther than point 1: keeps its value,
//   - point 3 (a node) at 0: of 2 and 3 only 2 is a frozen unknown (already
//     kept); point 3 is rounded,
//   - point 4 at 0.375: point 3 is no frozen coordinate, so point 4 keeps
//     its value; point 5 at 0.375 is as far as point 4 and is rounded,
//   - point 18 at 16777217 (not a float value), copy 12 at 8, point 6 at 4,
//   - copy 30: neither unknown nor eliminated - no entry.
//  Hole: its only kept point 14 is the neighbour of itself: u at 3.75 is a
//  frozen unknown and keeps its value, v (a node) at 3.75 is rounded to 4.
void borderRoundingSelection()
{
    parts::Part part;
    flatGrid(part.mesh);
    std::vector<int> corners = part.mesh.topology().cornerVertex();
    corners.erase(corners.begin() + 3 * 14, corners.begin() + 3 * 16);
    part.mesh.editTopology().swapCorners(corners, FaceMode::Triangles);
    part.mesh.topologyChanged();
    const BoundaryLoops loops(part.mesh);
    EXPECT(loops.count(), 2);
    EXPECT(v(loops.points().row(0)),
           (std::vector<int>{0, 1, 2, 3, 4, 5, 11, 17, 23, 29, 28, 27, 26, 25, 24, 18, 12, 6}));
    EXPECT(v(loops.points().row(1)), (std::vector<int>{9, 8, 14, 15}));

    const Adjacency& a = part.mesh.adjacency();
    field::CutGraph cuts;
    cuts.cutEdges = edgesBetween(a, {{12, 13}});
    cuts.cutEdgeBits = BitArray(static_cast<std::size_t>(a.edgeCount()));
    cuts.cutEdgeBits.set(static_cast<std::size_t>(cuts.cutEdges[0]));
    cuts.copy = field::CutCopy::cut(part.mesh, cuts.cutEdges, cuts.cutEdgeBits);
    part.cutGraph = std::move(cuts);
    freezeBorders(part, {0, 1, 2, 3, 4, 5, 18, 12, 6, 14});
    const param::SystemContext context(part, false);
    EXPECT(context.copy().pointCount(), 31);
    EXPECT(context.copies(12), (std::vector<int>{12, 30}));
    EXPECT(context.copy().adjacency().edgeBetween(30, 18) != -1 && context.copy().adjacency().edgeBetween(12, 6) != -1,
           true);

    param::System system;
    system.copyPointCount = 31;
    system.variableCount = 62;
    system.flags.assign(62, param::variable_flag::frozenBorder);
    for (const int variable : {36, 7, 29}) {
        system.flags[static_cast<std::size_t>(variable)] = param::variable_flag::node;
    }
    system.constraints = param::ConstraintSet(62);
    system.constraints.add(0, {});
    system.constraints.add(1, {});
    system.constraints.closeAnchors();
    system.constraints.add(60, {{24, 1.0}, {2, -0.125}});
    std::vector<int> integers;
    for (int variable = 0; variable < 62; ++variable) {
        const bool unknown = variable != 0 && variable != 1 && variable != 60 && variable != 61;
        system.variableToUnknown.push_back(unknown ? variable : -1);
        system.unknownToVariable.push_back(variable);
        if (unknown && variable != 8) {
            integers.push_back(variable);
        }
    }
    std::vector<double> x(62, 0.0);
    const std::initializer_list<std::pair<int, double>> values = {
        {2, 2.25}, {4, 2.375}, {6, 1.625}, {8, 4.5},   {10, -0.5}, {36, 7.0},         {24, 7.25}, {12, 0.25},
        {3, 0.125}, {5, -0.25}, {7, 0.0},  {9, 0.375}, {11, 0.375}, {37, 16777217.0}, {25, 8.0},  {13, 4.0},
        {28, 3.75}, {29, 3.75}};
    for (const auto& [variable, value] : values) {
        x[static_cast<std::size_t>(variable)] = value;
    }

    QuietProgress q;
    solver::CompactMatrix matrix;
    std::vector<double> rhs;
    NoAnswer callback;
    solver::ConstrainedDriver driver(matrix, rhs, callback, q.progress);
    driver.setIntegerUnknowns(integers);
    solver::DriverRequest request;
    request.x = &x;
    RecordedBorders borders;
    (void)borders.select(param::SolveState{context, system, driver, q.progress}, request);
    EXPECT(borders.outcome == param::RoundingStep::Outcome::Delivered, true);
    EXPECT(borders.picks, (std::vector<std::pair<int, double>>{
                              {2, 2.0},   {4, 2.375}, {6, 2.0}, {10, -1.0},  {36, 7.0}, {24, 7.0},
                              {12, 0.25}, {3, 0.125}, {5, -0.25}, {7, 0.0},  {9, 0.375}, {11, 0.0},
                              {37, 16777217.0}, {25, 8.0}, {13, 4.0}, {28, 3.75}, {29, 4.0}}));
    // Position 0 of the integer list is unknown 2, not unknown 0.
    EXPECT(request.selection.empty() ? -1 : request.selection.front().position, 0);
}

// The flat grid as a part with frozen borders and its system; `kept` are
// the border points that stay. Every coordinate of a border point has to
// become whole, the anchor (point 0) excepted: 34 unknowns.
parts::Part frozenFlatPart(std::initializer_list<int> kept)
{
    parts::Part part = flatPart({}, 0, {}, false);
    freezeBorders(part, kept);
    QuietProgress q;
    param::SystemBuilder::compute(part, param::SystemBuilder::Options{}, q.progress);
    return part;
}

// Where the rounding along frozen borders stands in the rounds, on the flat
// grid (the perfect grid u = 2x, v = 4y solves it). The regions are stood in
// for by a step that delivers the first position (u of point 1) at once.
//  - kept points 0 to 5 (the lower side): reweighting, regions, borders, the
//    rest - four rounds. The borders select 9 entries: u of the points 2 to
//    5, exactly 4, 6, 8, 10 (all different, none in conflict; u of point 1
//    is fixed already and left out), then v of the points 1 to 5 near 0,
//  - no kept point: the borders select nothing and everything left is
//    rounded in the same round - three rounds, every coordinate that had to
//    become whole is whole.
void borderRoundingRounds()
{
    QuietProgress q;
    {
        parts::Part part = frozenFlatPart({0, 1, 2, 3, 4, 5});
        EXPECT(part.system.has_value() && part.system->integerUnknowns.size() == 34, true);
        FirstPosition regions;
        RecordedBorders borders;
        param::SolutionBuilder::compute(part, param::SolutionBuilder::Options{}, q.progress,
                                        param::SolutionBuilder::Steps{regions, borders});
        EXPECT(part.solution.has_value(), true);
        if (!part.solution) {
            return;
        }
        EXPECT(regions.calls, 1);
        EXPECT(borders.calls, 1);
        EXPECT(borders.alreadyFixed, 0);
        EXPECT(part.solution->rounds, 4);
        std::vector<int> variables;
        std::vector<double> u;
        double worstV = 0.0;
        for (const auto& [variable, value] : borders.picks) {
            variables.push_back(variable);
            if (variable % 2 == 0) {
                u.push_back(value);
            } else {
                worstV = std::max(worstV, std::fabs(value));
            }
        }
        EXPECT(variables, (std::vector<int>{4, 6, 8, 10, 3, 5, 7, 9, 11}));
        EXPECT(u, (std::vector<double>{4.0, 6.0, 8.0, 10.0}));
        EXPECT(worstV < 1e-6, true);
        EXPECT(gridMisfit(part) < 1e-4, true);
    }
    {
        parts::Part part = frozenFlatPart({});
        EXPECT(part.system.has_value() && part.system->integerUnknowns.size() == 34, true);
        FirstPosition regions;
        RecordedBorders borders;
        param::SolutionBuilder::compute(part, param::SolutionBuilder::Options{}, q.progress,
                                        param::SolutionBuilder::Steps{regions, borders});
        EXPECT(part.solution.has_value(), true);
        if (!part.solution) {
            return;
        }
        EXPECT(borders.calls, 1);
        EXPECT(borders.outcome == param::RoundingStep::Outcome::Failed, true);
        EXPECT(borders.picks.empty(), true);
        EXPECT(part.solution->rounds, 3);
        EXPECT(notWhole(part), 0);
        EXPECT(gridMisfit(part) < 1e-4, true);
    }
}

// The round of the borders behind the rounding by regions of the core: the
// flat grid with its border as crease curves, frozen, the lower side kept
// (points 0 to 5). The curves tie u along the left and the right side and v
// along the lower and the upper side as in systemOnDisc; of the coordinates
// of the border points 16 are left as unknowns, and all of them have to
// become whole: u of the points 1 to 4 (variables 2 to 8), v of the points
// 6, 12 and 18 (13, 25, 37) and of the points 11, 17 and 23 (23, 35, 47), u
// of the points 25 to 28 (50 to 56) and the two free values of the curves,
// u of the right side (46) and v of the upper side (57).
//  - reweighting,
//  - regions: the two free values of the curves, as without frozen borders,
//  - borders: u of the kept points 1 to 4 lies at 2, 4, 6 and 8, no two
//    neighbours land on the same number, and each is selected with its whole
//    number. u of the points 0 and 5 and v of all six follow from conditions
//    and have no position to select,
//  - the ten unknowns left are rounded at once: four rounds.
void bordersBehindRegions()
{
    using Outcome = param::RoundingStep::Outcome;
    QuietProgress q;
    parts::Part part = flatPart({}, 0, {}, true);
    freezeBorders(part, {0, 1, 2, 3, 4, 5});
    param::SystemBuilder::compute(part, param::SystemBuilder::Options{}, q.progress);
    EXPECT(part.system.has_value(), true);
    if (!part.system) {
        return;
    }
    std::vector<int> whole;
    for (const int unknown : part.system->integerUnknowns) {
        whole.push_back(part.system->unknownToVariable[static_cast<std::size_t>(unknown)]);
    }
    EXPECT(whole, (std::vector<int>{2, 4, 6, 8, 13, 23, 25, 35, 37, 46, 47, 50, 52, 54, 56, 57}));

    RecordedRegions regions;
    RecordedBorders borders;
    param::SolutionBuilder::compute(part, param::SolutionBuilder::Options{}, q.progress,
                                    param::SolutionBuilder::Steps{regions, borders});
    EXPECT(part.solution.has_value(), true);
    if (!part.solution) {
        return;
    }
    EXPECT(regions.outcomes == (std::vector<Outcome>{Outcome::Delivered}), true);
    EXPECT(regions.selected, (std::vector<std::size_t>{2}));
    EXPECT(borders.calls, 1);
    EXPECT(borders.outcome == Outcome::Delivered, true);
    EXPECT(borders.alreadyFixed, 0);
    EXPECT(borders.picks, (std::vector<std::pair<int, double>>{{2, 2.0}, {4, 4.0}, {6, 6.0}, {8, 8.0}}));
    EXPECT(part.solution->rounds, 4);
    EXPECT(part.solution->integerUnknowns, 16);
    EXPECT(notWhole(part), 0);
    EXPECT(gridMisfit(part) < 1e-4, true);
}

// Entries of the rounding along frozen borders that the other checks do not
// meet, with system, marks and solution set by hand. The mesh is the flat
// grid without the squares 9 10 16 15 and 14 15 21 20: the two holes touch
// in point 15, which has four border edges, so the hole is walked as two
// loops that both hold point 15. No cut: copy point = point, variable 2p is
// u and 2p + 1 is v of point p, unknown = variable where there is one.
// Kept are the points 1 to 5 of the outer loop and 9 and 15.
//  - outer loop, u: point 1 at 2.25; point 2 is eliminated, its condition
//    names u of point 1 and a variable without unknown (u of point 25, at
//    100 in the solution vector), which adds nothing: both land on 2 and
//    point 1 keeps its value. Point 3 stands for a copy of the mirrored
//    half: no unknown, not a frozen coordinate, eliminated onto u of point 4
//    (3.75): both land on 4, so point 4 keeps its value. Point 5 has neither
//    unknown nor condition and gives no entry. v of the points 1 to 5 gives
//    no entry either,
//  - first hole loop (15 14 20 21): point 15 is its only kept point and its
//    own neighbour, so it keeps u at 6.125 and v at 8.375,
//  - second hole loop (15 16 10 9): u of point 15 at 6.125 and of point 9
//    at 5.75 both land on 6, point 9 is farther away and keeps its value;
//    point 15 would be rounded to 6 here and v of it to 8 (no conflict with
//    point 9 at 4), but both its positions are selected already and are
//    left out. Point 15 is selected once per coordinate.
void borderRoundingSharedPoint()
{
    parts::Part part;
    flatGrid(part.mesh);
    std::vector<int> corners = part.mesh.topology().cornerVertex();
    corners.erase(corners.begin() + 3 * 24, corners.begin() + 3 * 26);
    corners.erase(corners.begin() + 3 * 16, corners.begin() + 3 * 18);
    part.mesh.editTopology().swapCorners(corners, FaceMode::Triangles);
    part.mesh.topologyChanged();
    const BoundaryLoops loops(part.mesh);
    EXPECT(loops.count(), 3);
    if (loops.count() != 3) {
        return;
    }
    EXPECT(v(loops.points().row(0)),
           (std::vector<int>{0, 1, 2, 3, 4, 5, 11, 17, 23, 29, 28, 27, 26, 25, 24, 18, 12, 6}));
    EXPECT(v(loops.points().row(1)), (std::vector<int>{15, 14, 20, 21}));
    EXPECT(v(loops.points().row(2)), (std::vector<int>{15, 16, 10, 9}));

    const Adjacency& a = part.mesh.adjacency();
    field::CutGraph cuts;
    cuts.cutEdgeBits = BitArray(static_cast<std::size_t>(a.edgeCount()));
    cuts.copy = field::CutCopy::cut(part.mesh, cuts.cutEdges, cuts.cutEdgeBits);
    part.cutGraph = std::move(cuts);
    freezeBorders(part, {1, 2, 3, 4, 5, 9, 15});
    const param::SystemContext context(part, false);

    param::System system;
    system.copyPointCount = 30;
    system.variableCount = 60;
    system.flags.assign(60, param::variable_flag::continuous);
    for (const int variable : {2, 4, 8, 18, 19, 30, 31}) {
        system.flags[static_cast<std::size_t>(variable)] = param::variable_flag::frozenBorder;
    }
    system.constraints = param::ConstraintSet(60);
    system.constraints.add(0, {});
    system.constraints.add(1, {});
    system.constraints.closeAnchors();
    system.constraints.add(4, {{2, 1.0}, {50, 1.0}});
    system.constraints.add(6, {{8, 1.0}});
    const std::vector<int> without{0, 1, 3, 4, 5, 6, 7, 9, 10, 11, 50};
    for (int variable = 0; variable < 60; ++variable) {
        const bool unknown = std::find(without.begin(), without.end(), variable) == without.end();
        system.variableToUnknown.push_back(unknown ? variable : -1);
        system.unknownToVariable.push_back(variable);
    }
    std::vector<double> x(60, 0.0);
    const std::initializer_list<std::pair<int, double>> values = {{2, 2.25},  {8, 3.75},  {18, 5.75},
                                                                  {19, 4.0},  {30, 6.125}, {31, 8.375},
                                                                  {50, 100.0}};
    for (const auto& [variable, value] : values) {
        x[static_cast<std::size_t>(variable)] = value;
    }

    QuietProgress q;
    solver::CompactMatrix matrix;
    std::vector<double> rhs;
    NoAnswer callback;
    solver::ConstrainedDriver driver(matrix, rhs, callback, q.progress);
    driver.setIntegerUnknowns({2, 8, 18, 19, 30, 31});
    solver::DriverRequest request;
    request.x = &x;
    RecordedBorders borders;
    (void)borders.select(param::SolveState{context, system, driver, q.progress}, request);
    EXPECT(borders.outcome == param::RoundingStep::Outcome::Delivered, true);
    EXPECT(borders.picks, (std::vector<std::pair<int, double>>{
                              {2, 2.25}, {8, 3.75}, {30, 6.125}, {31, 8.375}, {18, 5.75}, {19, 4.0}}));
}

// The rounding along frozen borders on a part with border curves, in the
// situation of the measurement on open meshes: the flat grid with its border
// as crease curves and a first direction of size 0.48 (u = x / 0.48), frozen,
// the right side kept (points 5, 11, 17, 23, 29). u along the right side is
// tied by the iso lines to its reference, u of point 23 (variable 46); it
// is a frozen coordinate (kind 6, given before the iso lines and kept by
// them) and the only entry of the side with an unknown. All entries of the
// side land on 10 and neighbours are joined, so the reference is marked and
// keeps its value. v of the points 11, 17 and 23 (variables 23, 35, 47) lie
// at 4, 8 and 12.
//  - regions stood in for by a step that delivers the first position (u of
//    point 1): the reference is free and is selected with its value as it
//    is, near 10.38; the whole right side stays there, off the grid line,
//  - the regions of the core: they deliver the reference, rounded to 10,
//    before; the borders leave it out and the right side lies on 10.
void bordersWithCurves()
{
    using Outcome = param::RoundingStep::Outcome;
    QuietProgress q;
    for (const bool coreRegions : {false, true}) {
        parts::Part part = flatPart({}, 0, {}, true);
        part.sizeField->sizeU.assign(part.sizeField->sizeU.size(), 0.48f);
        freezeBorders(part, {5, 11, 17, 23, 29});
        param::SystemBuilder::compute(part, param::SystemBuilder::Options{}, q.progress);
        EXPECT(part.system.has_value(), true);
        if (!part.system) {
            return;
        }
        const param::System& system = *part.system;
        EXPECT(system.flags[46], param::variable_flag::frozenBorder);
        EXPECT(system.variableToUnknown[46] >= 0, true);
        for (const int variable : {10, 22, 34, 58}) {
            EXPECT(system.variableToUnknown[static_cast<std::size_t>(variable)], -1);
        }
        RecordedRegions coreStep;
        FirstPosition standIn;
        param::RoundingStep& regions = coreRegions ? static_cast<param::RoundingStep&>(coreStep) : standIn;
        RecordedBorders borders;
        param::SolutionBuilder::compute(part, param::SolutionBuilder::Options{}, q.progress,
                                        param::SolutionBuilder::Steps{regions, borders});
        EXPECT(part.solution.has_value(), true);
        if (!part.solution) {
            return;
        }
        EXPECT(part.solution->rounds, 4);
        EXPECT(borders.calls, 1);
        EXPECT(borders.alreadyFixed, 0);
        const param::UvOutput::Coordinates& uv = part.solution->uv;
        const bool sideEven = uv[5][0] == uv[11][0] && uv[5][0] == uv[17][0] && uv[5][0] == uv[23][0] &&
                              uv[5][0] == uv[29][0];
        EXPECT(sideEven, true);
        const std::vector<std::pair<int, double>> crossing{{23, 4.0}, {35, 8.0}, {47, 12.0}};
        if (coreRegions) {
            EXPECT(coreStep.outcomes == (std::vector<Outcome>{Outcome::Delivered}), true);
            EXPECT(borders.picks, crossing);
            EXPECT(uv[23][0], 10.0f);
            continue;
        }
        EXPECT(standIn.calls, 1);
        EXPECT(borders.picks.size(), std::size_t{4});
        if (borders.picks.size() != 4) {
            continue;
        }
        const auto [reference, value] = borders.picks.front();
        EXPECT(reference, 46);
        EXPECT(value > 10.3 && value < 10.4, true);
        EXPECT(uv[23][0], static_cast<float>(value));
        EXPECT((std::vector<std::pair<int, double>>(borders.picks.begin() + 1, borders.picks.end())), crossing);
    }
}

// Polygon splitting: points only, the polygon is 0 .. n-1.
Mesh ngonPoints(std::vector<Vec3f> points)
{
    Mesh m;
    m.swapPoints(points);
    return m;
}

std::vector<int> ngonRange(int n)
{
    std::vector<int> out(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        out[static_cast<std::size_t>(i)] = i;
    }
    return out;
}

using Polygons = std::vector<std::vector<int>>;

Polygons ngonPolygons(const extract::NgonPolygonList& list)
{
    Polygons out;
    for (int p = 0; p < list.polygonCount(); ++p) {
        out.push_back(v(list.polygon(p)));
    }
    return out;
}

std::vector<int> ngonFlags(std::vector<Vec3f> points)
{
    const int n = static_cast<int>(points.size());
    Mesh m = ngonPoints(std::move(points));
    extract::NgonNewPoints fresh(m);
    extract::NgonPolygonList list(ngonRange(n));
    extract::NgonCornerData data;
    extract::CornerClassifier::classify(list, fresh, data);
    return data.flag;
}

Polygons ngonSplit(std::vector<Vec3f> points, int* added = nullptr)
{
    const int n = static_cast<int>(points.size());
    Mesh m = ngonPoints(std::move(points));
    extract::NgonNewPoints fresh(m);
    extract::NgonPolygonList list(ngonRange(n));
    extract::PolygonSplitter::split(list, fresh);
    if (added != nullptr) {
        *added = fresh.count();
    }
    return ngonPolygons(list);
}

Polygons meshFaces(const Mesh& m)
{
    Polygons out;
    const Topology& t = m.topology();
    for (int f = 0; f < t.faceCount(); ++f) {
        out.emplace_back(t.cornerVertex().begin() + t.faceBegin(f), t.cornerVertex().begin() + t.faceEnd(f));
    }
    return out;
}

const std::vector<Vec3f> chamferedSquare{{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 2, 0}, {0, 2, 0}};
const std::vector<Vec3f> lHexagon{{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}, {1, 2, 0}, {0, 2, 0}};
const std::vector<Vec3f> sharpHexagon{{0, 0, 0}, {-2, 2, 0}, {-3, 1, 0}, {-2, 0, 0}, {-1.5f, 0, 0}, {-0.75f, 0, 0}};
const std::vector<Vec3f> triangleWithMiddles{{0, 0, 0},    {1, 0, 0}, {2, 0, 0},
                                             {1.5f, 0.9f, 0}, {1, 1.8f, 0}, {0.5f, 0.9f, 0}};
const std::vector<Vec3f> stripNonagon{{0, 1, 0},    {0, 0, 0},     {0.5f, 0, 0},  {2.1f, 0, 0},  {3, 0, 0},
                                      {3, 1, 0},    {2.25f, 1, 0}, {1.5f, 1, 0},  {0.75f, 1, 0}};

void ngonCornerFlags()
{
    // Corner, straight, corner, corner, corner.
    EXPECT(ngonFlags({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {0, 1, 0}}), (std::vector<int>{1, 0, 1, 1, 1}));
    // One reflex corner.
    EXPECT(ngonFlags(lHexagon), (std::vector<int>{1, 1, 1, -1, 1, 1}));
    // Three corners and two turns of 45 degrees: the first 45 becomes the fourth
    // corner, the other one straight.
    EXPECT(ngonFlags(chamferedSquare), (std::vector<int>{1, 1, 1, 0, 1}));
    // Eight turns of 45 degrees: all straight.
    EXPECT(ngonFlags({{1, 0, 0}, {2, 0, 0}, {3, 1, 0}, {3, 2, 0}, {2, 3, 0}, {1, 3, 0}, {0, 2, 0}, {0, 1, 0}}),
           (std::vector<int>(8, 0)));
    // Turns of 135 and 45 degrees stay unclear (two corners, sum far from a full turn).
    EXPECT(ngonFlags(sharpHexagon), (std::vector<int>{3, 1, 1, 3, 0, 0}));
    // With four corners and unclear ones, the second look makes them corners.
    {
        Mesh m = ngonPoints(sharpHexagon);
        extract::NgonNewPoints fresh(m);
        extract::NgonPolygonList list(ngonRange(6));
        extract::NgonCornerData data;
        extract::CornerClassifier::classify(list, fresh, data);
        extract::FlagCounts k = extract::CornerClassifier::count(data, 0, 6);
        EXPECT(k.corners, 2);
        EXPECT(k.unclear, 2);
        EXPECT(k.lastStraight, 5);
        extract::CornerClassifier::reclassify(data, k, 0, 6);
        EXPECT(data.flag, (std::vector<int>{1, 1, 1, 1, 0, 0}));
        EXPECT(k.corners, 4);
        EXPECT(k.unclear, 0);
    }
}

void ngonRules()
{
    int added = -1;
    // Pentagon with four corners and one straight: chord from the straight to the corner opposite.
    EXPECT(ngonSplit(chamferedSquare, &added), (Polygons{{0, 1, 3, 4}, {1, 2, 3}}));
    EXPECT(added, 0);
    // Pentagon with a spike: the cheapest chord i -> i+3 (short, square to its straight end).
    EXPECT(ngonSplit({{0, 0, 0}, {2, 0, 0}, {4, 0, 0}, {2, 0.5f, 0}, {0, 0.5f, 0}}), (Polygons{{0, 1, 3, 4}, {1, 2, 3}}));
    // Reflex corner: chord from it.
    EXPECT(ngonSplit(lHexagon), (Polygons{{0, 3, 4, 5}, {0, 1, 2, 3}}));
    // Pentagon with a reflex corner: planar triangulation (a fan here, with a flat triangle).
    EXPECT(ngonSplit({{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {1, 1, 0}, {0, 2, 0}}),
           (Polygons{{0, 1, 2}, {0, 2, 3}, {0, 3, 4}}));
    // Sides of 3 and 5 edges between two corner pairs: three quads.
    EXPECT(ngonSplit({{0, 1, 0}, {0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1.5f, 1, 0}, {1, 1, 0}, {0.5f, 1, 0}}),
           (Polygons{{2, 3, 4, 5}, {5, 6, 7, 2}, {7, 0, 1, 2}}));
    // Cut-off quad at an unclear corner.
    EXPECT(ngonSplit({{0, 0, 0}, {4, 0, 0}, {4, 0.985f, 0}, {3, 0.985f, 0}, {1.706f, 0.985f, 0}, {0.766f, 0.643f, 0}}),
           (Polygons{{0, 3, 4, 5}, {0, 1, 2, 3}}));
    // Three corners between straights: one new point inside.
    {
        Mesh m = ngonPoints({{1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {0, 1, 0}, {0, 0.5f, 0}, {0, 0, 0}});
        extract::NgonNewPoints fresh(m);
        extract::NgonPolygonList list(ngonRange(6));
        extract::PolygonSplitter::split(list, fresh);
        EXPECT(ngonPolygons(list), (Polygons{{0, 1, 2, 6}, {2, 3, 4, 6}, {4, 5, 0, 6}}));
        EXPECT(fresh.positions(), (std::vector<Vec3f>{{1, 0.5f, 0}}));
        EXPECT(fresh.parents().entries, (std::vector<int>{2, 0, 1, 4, 3}));
        EXPECT(fresh.weights(), (std::vector<float>{1, 1, -1, 1, -1}));
    }
    // Corner, two straights, corner: two new points cut off the quads at both corners.
    {
        Mesh m = ngonPoints(sharpHexagon);
        extract::NgonNewPoints fresh(m);
        extract::NgonPolygonList list(ngonRange(6));
        extract::PolygonSplitter::split(list, fresh);
        EXPECT(ngonPolygons(list), (Polygons{{4, 5, 7, 6}, {0, 1, 7, 5}, {6, 7, 1, 2}, {3, 4, 6, 2}}));
        EXPECT(fresh.parents().starts, (std::vector<int>{0, 3, 6}));
        EXPECT(fresh.parents().entries, (std::vector<int>{4, 2, 3, 5, 1, 0}));
        EXPECT(fresh.weights(), (std::vector<float>{1, 0.6f, -0.6f, 1, 0.6f, -0.6f}));
    }
    // Single corner between straights: new point opposite the corner.
    {
        Mesh m = ngonPoints(triangleWithMiddles);
        extract::NgonNewPoints fresh(m);
        extract::NgonPolygonList list(ngonRange(6));
        extract::PolygonSplitter::split(list, fresh);
        EXPECT(ngonPolygons(list), (Polygons{{5, 0, 1, 6}, {1, 2, 3, 6}, {3, 4, 5, 6}}));
        EXPECT(fresh.parents().entries, (std::vector<int>{0, 5, 1}));
        EXPECT(fresh.weights(), (std::vector<float>{-1, 1, 1}));
    }
}

void ngonStripAndGrid()
{
    // Strip: the first two corner pairs give sides of 4 and 5 edges; the
    // chord with the closest arc shares joins position 3 and 6.
    {
        Mesh m = ngonPoints(stripNonagon);
        extract::NgonNewPoints fresh(m);
        extract::NgonPolygonList list(ngonRange(9));
        extract::NgonCornerData data;
        extract::CornerClassifier::classify(list, fresh, data);
        EXPECT(data.flag, (std::vector<int>{1, 1, 0, 0, 1, 1, 0, 0, 0}));
        extract::NgonSplitWork work{list, data, fresh};
        EXPECT(extract::StripGrid::strip(work, 0), true);
        EXPECT(ngonPolygons(list), (Polygons{{0, 1, 2, 3, 6, 7, 8}, {3, 4, 5, 6}}));
    }
    // Grid: a 2 x 2 outline gets one inner point in the middle (unchanged by the smoothing).
    {
        Mesh m = ngonPoints({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {2, 2, 0}, {1, 2, 0}, {0, 2, 0}, {0, 1, 0}});
        extract::NgonNewPoints fresh(m);
        extract::NgonPolygonList list(ngonRange(8));
        extract::PolygonSplitter::split(list, fresh);
        EXPECT(ngonPolygons(list), (Polygons{{0, 1, 8, 7}, {1, 2, 3, 8}, {7, 8, 5, 6}, {8, 3, 4, 5}}));
        EXPECT(fresh.positions(), (std::vector<Vec3f>{{1, 1, 0}}));
        EXPECT(fresh.parents().entries, (std::vector<int>{7, 3, 1, 5}));
        EXPECT(fresh.weights(), (std::vector<float>{0.25f, 0.25f, 0.25f, 0.25f}));
    }
}

void ngonNewPointParents()
{
    Mesh m = ngonPoints(grid(3));
    extract::NgonNewPoints fresh(m);
    const int x = fresh.add({0, 0, 0});
    fresh.addParent(0, 1.0f);
    fresh.addParent(1, 0.5f);
    fresh.endParents();
    const int y = fresh.add({0, 0, 0});
    fresh.addParent(x, 0.5f);
    fresh.addParent(2, 1.0f);
    fresh.endParents();
    EXPECT(x, 3);
    EXPECT(y, 4);
    EXPECT(fresh.parents().starts, (std::vector<int>{0, 2, 5}));
    EXPECT(fresh.parents().entries, (std::vector<int>{0, 1, 0, 1, 2}));
    EXPECT(fresh.weights(), (std::vector<float>{1.0f, 0.5f, 0.5f, 0.25f, 1.0f}));
    // A new round starts without parents and weights.
    fresh.reset();
    EXPECT(fresh.count(), 0);
    EXPECT(fresh.weights().empty(), true);
    EXPECT(fresh.parents().starts, (std::vector<int>{0}));
}

void ngonParity()
{
    // Row: triangle, quad, quad, triangle. Middle points on the three shared
    // edges; the end faces gain a corner, the quads are cut lengthwise.
    {
        Mesh m = polygonMesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}, {0, 1, 0}, {-1, 0.5f, 0}, {3, 0.5f, 0}},
                             {0, 3, 7, 11, 14}, {6, 0, 5, 0, 1, 4, 5, 1, 2, 3, 4, 2, 7, 3}, FaceMode::Polygons);
        extract::SplitRound round(m);
        extract::ParityBalance::run(round, BitArray(4, true));
        EXPECT(round.pendingCount(), 6);
        round.commit();
        EXPECT(meshFaces(m), (Polygons{{6, 0, 8, 5}, {8, 9, 4, 5}, {9, 8, 0, 1}, {9, 10, 3, 4}, {10, 9, 1, 2}, {10, 2, 7, 3}}));
        EXPECT(std::vector<Vec3f>(m.points().begin() + 8, m.points().end()),
               (std::vector<Vec3f>{{0, 0.5f, 0}, {1, 0.5f, 0}, {2, 0.5f, 0}}));
    }
    // Triangle left of a quad, triangle below it: the path turns at the shared
    // corner, a point inside joins both middles.
    {
        Mesh m = polygonMesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {-1, 0.5f, 0}, {0.5f, -1, 0}}, {0, 3, 7, 10},
                             {4, 0, 3, 0, 1, 2, 3, 0, 5, 1}, FaceMode::Polygons);
        extract::SplitRound round(m);
        extract::ParityBalance::run(round, BitArray(3, true));
        EXPECT(round.points().parents().starts, (std::vector<int>{0, 2, 4, 7}));
        EXPECT(round.points().parents().entries, (std::vector<int>{0, 3, 0, 1, 3, 1, 0}));
        EXPECT(round.points().weights(), (std::vector<float>{0.5f, 0.5f, 0.5f, 0.5f, 0.4f, 0.4f, 0.2f}));
        round.commit();
        EXPECT(meshFaces(m), (Polygons{{4, 0, 6, 3}, {0, 7, 8, 6}, {8, 7, 1, 2}, {8, 2, 3, 6}, {7, 0, 5, 1}}));
        EXPECT(std::vector<Vec3f>(m.points().begin() + 6, m.points().end()),
               (std::vector<Vec3f>{{0, 0.5f, 0}, {0.5f, 0, 0}, {0.4f, 0.4f, 0}}));
    }
}

void ngonDoubletsAndSmoothing()
{
    // A point between two edges inside a square split into two quads.
    Mesh m = polygonMesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5f, 0.5f, 0}}, {0, 4, 8},
                         {0, 1, 2, 4, 2, 3, 0, 4}, FaceMode::Polygons);
    const Adjacency& a = m.adjacency();
    std::vector<int> expected{a.edgeBetween(4, 0), a.edgeBetween(4, 2)};
    std::sort(expected.begin(), expected.end());
    EXPECT(extract::DoubletEdges::find(m), expected);
    // Dissolved once: one quad, the doublet point stays without corners.
    extract::NgonOptimizer::run(m, extract::NgonSettings::forPart(false, false));
    EXPECT(m.faceCount(), 1);
    EXPECT(m.topology().faceSize(0), 4);
    EXPECT(m.pointCount(), 5);
    EXPECT(m.adjacency().pointCorners(4).empty(), true);
    // Not a doublet when one face is no quad.
    Mesh five = polygonMesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5f, 0.5f, 0}, {0.5f, 0, 0}}, {0, 5, 9},
                            {0, 5, 1, 2, 4, 2, 3, 0, 4}, FaceMode::Polygons);
    EXPECT(extract::DoubletEdges::find(five).empty(), true);

    // Smoothing: point 4 inside a square of side 2 (four triangles); point 5 has no edges and stays.
    Mesh fan = polygonMesh({{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0.5f, 0.5f, 0}, {7, 7, 0}}, {0, 3, 6, 9, 12},
                           {0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4}, FaceMode::Polygons);
    extract::NewPointSmoother::run(fan, 4, 2, 0.5f);
    EXPECT(fan.points()[4], (Vec3f{0.875f, 0.875f, 0}));
    EXPECT(fan.points()[5], (Vec3f{7, 7, 0}));
}

void ngonWholeRun()
{
    auto check = [](const Mesh& m) {
        int largest = 0;
        for (int f = 0; f < m.faceCount(); ++f) {
            largest = std::max(largest, m.topology().faceSize(f));
        }
        EXPECT(largest <= 4, true);
        EXPECT(m.adjacency().nonManifoldEdges().empty(), true);
    };
    // Four separate polygons, without parity rounds, with smoothing.
    std::vector<Vec3f> points;
    std::vector<int> starts{0};
    std::vector<int> corners;
    float shift = 0.0f;
    for (const auto* polygon : {&triangleWithMiddles, &chamferedSquare, &stripNonagon, &lHexagon}) {
        const int first = static_cast<int>(points.size());
        for (const Vec3f& p : *polygon) {
            points.push_back({p.x + shift, p.y, p.z});
            corners.push_back(first + static_cast<int>(&p - polygon->data()));
        }
        starts.push_back(static_cast<int>(corners.size()));
        shift += 10.0f;
    }
    Mesh m = polygonMesh(points, starts, corners, FaceMode::Polygons);
    extract::NgonOptimizer::run(m, extract::NgonSettings::forPart(false, true));
    check(m);
    // Triangles and a quad with parity rounds: quads only.
    Mesh l = polygonMesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {-1, 0.5f, 0}, {0.5f, -1, 0}}, {0, 3, 7, 10},
                         {4, 0, 3, 0, 1, 2, 3, 0, 5, 1}, FaceMode::Polygons);
    extract::NgonOptimizer::run(l, extract::NgonSettings::forPart(true, false));
    check(l);
    EXPECT(l.faceCount(), 5);
    EXPECT(l.adjacency().counts().quads, 5);
}

// Two commits with parity rounds: the even octagon (points 6 .. 13) first
// gets its inner point 14 with four weights; then the parity round of the
// two triangles beside the quad adds the middles 15, 16 and the inner point
// 17 = 0.4 p3 + 0.4 p1 + 0.2 p0, and the next round, which looks only at the
// five new faces, finds nothing. A point map follows the weights of the
// second round (12), not equal shares (10).
void ngonTwoCommits()
{
    Mesh m = polygonMesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {-1, 0.5f, 0}, {0.5f, -1, 0}, {10, 0, 0},
                          {11, 0, 0}, {12, 0, 0}, {12, 1, 0}, {12, 2, 0}, {11, 2, 0}, {10, 2, 0}, {10, 1, 0}},
                         {0, 3, 7, 10, 18}, {4, 0, 3, 0, 1, 2, 3, 0, 5, 1, 6, 7, 8, 9, 10, 11, 12, 13},
                         FaceMode::Polygons);
    std::vector<float> values(14, 0.0f);
    values[1] = 10.0f;
    values[3] = 20.0f;
    m.maps().set(std::make_unique<TypedMap<float>>(std::string(map_names::density), ElementKind::Point, values));
    extract::NgonOptimizer::run(m, extract::NgonSettings::forPart(true, false));
    EXPECT(m.faceCount(), 9);
    EXPECT(m.adjacency().counts().quads, 9);
    EXPECT(m.pointCount(), 18);
    const auto* density = m.maps().findAs<TypedMap<float>>(map_names::density, ElementKind::Point);
    EXPECT(density != nullptr && density->size() == 18, true);
    if (density != nullptr && density->size() == 18) {
        EXPECT(std::fabs(density->values()[15] - 10.0f) < 1e-4f, true);
        EXPECT(std::fabs(density->values()[16] - 5.0f) < 1e-4f, true);
        EXPECT(std::fabs(density->values()[17] - 12.0f) < 1e-4f, true);
    }
}

} // namespace

namespace {

// A triangle mesh with grid coordinates per corner, for the trace.
struct GridMesh {
    std::vector<Vec3f> points;
    std::vector<int> corners;
    std::vector<std::array<float, 2>> cornerUv;

    void triangle(std::array<int, 3> p, std::array<std::array<float, 2>, 3> uv)
    {
        corners.insert(corners.end(), p.begin(), p.end());
        cornerUv.insert(cornerUv.end(), uv.begin(), uv.end());
    }
};

// The part of such a mesh: cut along the point pairs `cut`, with `jumpAt(e)`
// quarter turns left across cut edge e, the grid coordinates of the corners
// on the points of the cut copy, and the start points.
parts::Part tracedPart(const GridMesh& grid, std::initializer_list<std::pair<int, int>> cut,
                       const std::function<int(const Adjacency&, int)>& jumpAt, std::vector<int> singular,
                       std::vector<int> seeds)
{
    parts::Part part;
    std::vector<Vec3f> points = grid.points;
    std::vector<int> corners = grid.corners;
    part.mesh.swapPoints(points);
    part.mesh.editTopology().swapCorners(corners, FaceMode::Triangles);
    part.mesh.topologyChanged();
    const Adjacency& a = part.mesh.adjacency();
    const auto edgeCount = static_cast<std::size_t>(a.edgeCount());
    field::CutGraph cuts;
    cuts.cutEdgeBits = BitArray(edgeCount);
    cuts.residualJump.assign(edgeCount, 0);
    for (const int e : edgesBetween(a, cut)) {
        cuts.cutEdgeBits.set(static_cast<std::size_t>(e));
        cuts.residualJump[static_cast<std::size_t>(e)] = jumpAt(a, e);
    }
    cuts.cutEdgeBits.forEachSet([&](std::size_t e) { cuts.cutEdges.push_back(static_cast<int>(e)); });
    cuts.copy = field::CutCopy::cut(part.mesh, cuts.cutEdges, cuts.cutEdgeBits);
    param::Solution solution;
    solution.uv.assign(static_cast<std::size_t>(cuts.copy.pointCount()), {});
    for (int c = 0; c < part.mesh.cornerCount(); ++c) {
        solution.uv[static_cast<std::size_t>(cuts.copy.topology().vertex(c))] = grid.cornerUv[static_cast<std::size_t>(c)];
    }
    solution.singularities.points = std::move(singular);
    solution.singularities.index.assign(solution.singularities.points.size(), 3);
    solution.singularities.fullTurn.assign(solution.singularities.points.size(), 0);
    solution.seeds = std::move(seeds);
    part.cutGraph = std::move(cuts);
    part.solution = std::move(solution);
    return part;
}

// What the trace found, in grid points: `keyOf` names the grid point a node
// position lies on (-1 for none).
struct TracedGrid {
    std::array<int, 3> kinds{};
    int nodes = 0;
    int misplaced = 0;
    // Every face as its sorted grid points (four, else marked with -2).
    std::vector<std::array<int, 4>> faces;
};

TracedGrid traceGrid(parts::Part& part, const std::function<int(Vec3f)>& keyOf)
{
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    extract::LoopPolygons closer;
    extract::IsolineTracer::trace(part, closer, progress);
    TracedGrid out;
    if (!part.trace) {
        out.misplaced = -1;
        return out;
    }
    const extract::Trace& trace = *part.trace;
    out.nodes = trace.nodes.count();
    std::vector<int> keyOfId(static_cast<std::size_t>(trace.nodes.idCount()), -1);
    for (int n = 0; n < trace.nodes.count(); ++n) {
        const extract::TraceNode& node = trace.nodes[n];
        ++out.kinds[static_cast<std::size_t>(node.kind)];
        const int key = keyOf(node.position);
        out.misplaced += key < 0 ? 1 : 0;
        if (node.id >= 0) {
            keyOfId[static_cast<std::size_t>(node.id)] = key;
        }
    }
    for (int r = 0; r < trace.polygons.rows(); ++r) {
        const auto polygon = trace.polygons.row(r);
        std::array<int, 4> face{-2, -2, -2, -2};
        if (polygon.size() == 4) {
            for (std::size_t i = 0; i < 4; ++i) {
                face[i] = keyOfId[static_cast<std::size_t>(polygon[i])];
            }
            std::sort(face.begin(), face.end());
        }
        out.faces.push_back(face);
    }
    std::sort(out.faces.begin(), out.faces.end());
    return out;
}

std::vector<std::array<int, 4>> sortedCells(std::vector<std::array<int, 4>> cells)
{
    for (auto& cell : cells) {
        std::sort(cell.begin(), cell.end());
    }
    std::sort(cells.begin(), cells.end());
    return cells;
}

// The square [0, 3]^2 from a lattice of step 0.75 (cells split along the
// rising diagonal), (u, v) = (x, y), one start point in the corner (0, 0).
// 16 nodes: the start node on its point; the other three corners (on
// points of the mesh, which a new node only keeps until it is put onto the
// nearest edge of its triangle), the eight border nodes and (1, 1), (2, 2)
// on edges (the diagonals run through them); (1, 2) and (2, 1) inside
// triangles. The nine unit cells are the faces, each once.
void traceSquare()
{
    GridMesh grid;
    for (int j = 0; j < 5; ++j) {
        for (int i = 0; i < 5; ++i) {
            grid.points.push_back({0.75f * static_cast<float>(i), 0.75f * static_cast<float>(j), 0.0f});
        }
    }
    auto uv = [&](int p) { return std::array<float, 2>{grid.points[static_cast<std::size_t>(p)].x, grid.points[static_cast<std::size_t>(p)].y}; };
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            const int a = 5 * j + i;
            grid.triangle({a, a + 1, a + 6}, {uv(a), uv(a + 1), uv(a + 6)});
            grid.triangle({a, a + 6, a + 5}, {uv(a), uv(a + 6), uv(a + 5)});
        }
    }
    parts::Part part = tracedPart(grid, {}, [](const Adjacency&, int) { return 0; }, {}, {0});
    const TracedGrid traced = traceGrid(part, [](Vec3f p) {
        const float x = std::round(p.x);
        const float y = std::round(p.y);
        const bool on = std::fabs(p.x - x) < 1e-5f && std::fabs(p.y - y) < 1e-5f && p.z == 0.0f;
        return on ? static_cast<int>(x) + 4 * static_cast<int>(y) : -1;
    });
    EXPECT(traced.nodes, 16);
    EXPECT(traced.misplaced, 0);
    EXPECT(traced.kinds, (std::array<int, 3>{2, 1, 13}));
    std::vector<std::array<int, 4>> cells;
    for (int y = 0; y < 3; ++y) {
        for (int x = 0; x < 3; ++x) {
            const int k = x + 4 * y;
            cells.push_back({k, k + 1, k + 5, k + 4});
        }
    }
    EXPECT(traced.faces, sortedCells(cells));
}

// Three rhombi of 120 degrees around a centre point, each a lattice of step
// 2/3 in its own coordinates (s, t) from 0 to 2, (u, v) = (s, t), (-t, s),
// (-s, -t) in the sectors 0, 1, 2. The spoke between sector 2 and sector 0
// is cut: across it the grid turns by a quarter (one quarter turn left from
// a triangle of sector 0 to one of sector 2). The centre is the singular
// point. 19 nodes (centre, two per spoke, four per sector), all but the
// start node on edges, the 12 cells of the sectors are the faces, three of
// them at the centre.
void traceStar()
{
    const std::array<Vec3f, 3> spoke{Vec3f{0.0f, 1.0f, 0.0f}, Vec3f{-0.8660254f, -0.5f, 0.0f},
                                     Vec3f{0.8660254f, -0.5f, 0.0f}};
    auto at = [&](int i, float s, float t) {
        const Vec3f a = spoke[static_cast<std::size_t>(i)];
        const Vec3f b = spoke[static_cast<std::size_t>((i + 1) % 3)];
        return Vec3f{s * a.x + t * b.x, s * a.y + t * b.y, s * a.z + t * b.z};
    };
    // Point of lattice place (si, ti) of sector i.
    auto point = [](int i, int si, int ti) {
        if (si == 0 && ti == 0) {
            return 0;
        }
        if (ti == 0) {
            return 1 + 3 * i + si - 1;
        }
        if (si == 0) {
            return 1 + 3 * ((i + 1) % 3) + ti - 1;
        }
        return 10 + 9 * i + 3 * (si - 1) + ti - 1;
    };
    auto step = [](int k) { return static_cast<float>(2 * k) / 3.0f; };
    GridMesh grid;
    grid.points.resize(37);
    for (int i = 0; i < 3; ++i) {
        for (int si = 0; si < 4; ++si) {
            for (int ti = 0; ti < 4; ++ti) {
                grid.points[static_cast<std::size_t>(point(i, si, ti))] = at(i, step(si), step(ti));
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        auto uv = [&](int si, int ti) {
            const float s = step(si);
            const float t = step(ti);
            return i == 0 ? std::array<float, 2>{s, t}
                          : (i == 1 ? std::array<float, 2>{0.0f - t, s} : std::array<float, 2>{0.0f - s, 0.0f - t});
        };
        for (int si = 0; si < 3; ++si) {
            for (int ti = 0; ti < 3; ++ti) {
                grid.triangle({point(i, si, ti), point(i, si + 1, ti), point(i, si + 1, ti + 1)},
                              {uv(si, ti), uv(si + 1, ti), uv(si + 1, ti + 1)});
                grid.triangle({point(i, si, ti), point(i, si + 1, ti + 1), point(i, si, ti + 1)},
                              {uv(si, ti), uv(si + 1, ti + 1), uv(si, ti + 1)});
            }
        }
    }
    // Sector 0 owns the triangles 0 .. 17.
    auto jump = [](const Adjacency& a, int e) { return a.edgeCorners(e)[0] / 3 < 18 ? 1 : -1; };
    parts::Part part = tracedPart(grid, {{0, 1}, {1, 2}, {2, 3}}, jump, {0}, {});

    // Grid points: centre 0, spoke i at 1 and 2: 1 + 2i + k - 1, sector i
    // at (s, t) in {1, 2}^2: 7 + 4i + 2(s - 1) + t - 1.
    std::vector<Vec3f> gridPoints(19);
    for (int i = 0; i < 3; ++i) {
        for (int s = 0; s < 3; ++s) {
            for (int t = 0; t < 3; ++t) {
                const int key = s == 0 && t == 0 ? 0
                                : t == 0         ? 1 + 2 * i + s - 1
                                : s == 0         ? 1 + 2 * ((i + 1) % 3) + t - 1
                                                 : 7 + 4 * i + 2 * (s - 1) + t - 1;
                gridPoints[static_cast<std::size_t>(key)] = at(i, static_cast<float>(s), static_cast<float>(t));
            }
        }
    }
    const TracedGrid traced = traceGrid(part, [&](Vec3f p) {
        for (std::size_t k = 0; k < gridPoints.size(); ++k) {
            const Vec3f d = p - gridPoints[k];
            if (dot(d, d) < 1e-10f) {
                return static_cast<int>(k);
            }
        }
        return -1;
    });
    EXPECT(traced.nodes, 19);
    EXPECT(traced.misplaced, 0);
    EXPECT(traced.kinds, (std::array<int, 3>{0, 1, 18}));
    std::vector<std::array<int, 4>> cells;
    auto key = [](int i, int s, int t) {
        return s == 0 && t == 0 ? 0
               : t == 0         ? 1 + 2 * i + s - 1
               : s == 0         ? 1 + 2 * ((i + 1) % 3) + t - 1
                                : 7 + 4 * i + 2 * (s - 1) + t - 1;
    };
    for (int i = 0; i < 3; ++i) {
        for (int s = 0; s < 2; ++s) {
            for (int t = 0; t < 2; ++t) {
                cells.push_back({key(i, s, t), key(i, s + 1, t), key(i, s + 1, t + 1), key(i, s, t + 1)});
            }
        }
    }
    EXPECT(traced.faces, sortedCells(cells));
    EXPECT(std::count_if(traced.faces.begin(), traced.faces.end(), [](const auto& f) { return f[0] == 0; }), 3);
}

// An open cylinder of circumference 4 and height 2 from a lattice of six
// columns and three rows (step 2/3 in u and v), cut along the column at
// u = 0 (u = 4 on the other side, no turn), one start point on the lower
// border at the seam. 12 nodes, all but the start node on edges (u = 1 and
// 3 on the chords halfway between two columns), the 8 cells are the faces,
// two of them across the seam.
void traceCylinder()
{
    const float radius = 4.0f / (2.0f * std::numbers::pi_v<float>);
    auto at = [&](float u, float v) {
        const float angle = 2.0f * std::numbers::pi_v<float> * u / 4.0f;
        return Vec3f{radius * std::cos(angle), radius * std::sin(angle), v};
    };
    auto step = [](int k) { return static_cast<float>(2 * k) / 3.0f; };
    GridMesh grid;
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 6; ++col) {
            grid.points.push_back(at(step(col), step(row)));
        }
    }
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 6; ++col) {
            const int next = (col + 1) % 6;
            const float u0 = step(col);
            const float u1 = step(col + 1);
            const std::array<float, 2> a{u0, step(row)};
            const std::array<float, 2> b{u1, step(row)};
            const std::array<float, 2> c{u1, step(row + 1)};
            const std::array<float, 2> d{u0, step(row + 1)};
            const int pa = 6 * row + col;
            const int pb = 6 * row + next;
            const int pc = 6 * (row + 1) + next;
            const int pd = 6 * (row + 1) + col;
            grid.triangle({pa, pb, pc}, {a, b, c});
            grid.triangle({pa, pc, pd}, {a, c, d});
        }
    }
    parts::Part part = tracedPart(grid, {{0, 6}, {6, 12}, {12, 18}}, [](const Adjacency&, int) { return 0; }, {}, {0});
    // Grid point (u, v) on the mesh: on the chord between two columns of the
    // lattice (u = 1 and 3 lie halfway between two of them).
    auto onMesh = [&](int u, int v) {
        const int col = 3 * u / 2;
        const Vec3f a = at(step(col), static_cast<float>(v));
        const Vec3f b = at(step(col + 1), static_cast<float>(v));
        return u % 2 == 0 ? a : Vec3f{(a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f, a.z};
    };
    const TracedGrid traced = traceGrid(part, [&](Vec3f p) {
        for (int v = 0; v < 3; ++v) {
            for (int u = 0; u < 4; ++u) {
                const Vec3f d = p - onMesh(u, v);
                if (dot(d, d) < 1e-10f) {
                    return u + 4 * v;
                }
            }
        }
        return -1;
    });
    EXPECT(traced.nodes, 12);
    EXPECT(traced.misplaced, 0);
    EXPECT(traced.kinds, (std::array<int, 3>{0, 1, 11}));
    std::vector<std::array<int, 4>> cells;
    for (int v = 0; v < 2; ++v) {
        for (int u = 0; u < 4; ++u) {
            const int next = (u + 1) % 4;
            cells.push_back({u + 4 * v, next + 4 * v, next + 4 * (v + 1), u + 4 * (v + 1)});
        }
    }
    EXPECT(traced.faces, sortedCells(cells));
    EXPECT(std::count_if(traced.faces.begin(), traced.faces.end(),
                         [](const auto& f) { return std::count_if(f.begin(), f.end(), [](int k) { return k % 4 == 0; }) == 2 &&
                                                    std::count_if(f.begin(), f.end(), [](int k) { return k % 4 == 3; }) == 2; }),
           2);
}

// An open fan of five triangles k = (0, k + 1, k + 2) around point 0 (the
// outer points on a half circle), cut along the spokes 0-2 and 0-4 with one
// quarter turn left from triangle 1 to 0 and from triangle 2 to 3.
// Triangles 0 and 2 meet only in point 0. From triangle 2 the fan first runs
// over the spoke 0-4 into the border; the fan the other way round starts
// again from the given direction and reaches triangle 0 over the spoke 0-2:
// direction 0 there is direction 3 here.
void traceFanOtherWay()
{
    GridMesh grid;
    grid.points.push_back({0.0f, 0.0f, 0.0f});
    for (int k = 0; k < 6; ++k) {
        const float angle = std::numbers::pi_v<float> * static_cast<float>(k) / 5.0f;
        grid.points.push_back({std::cos(angle), std::sin(angle), 0.0f});
    }
    auto uv = [&](int p) { return std::array<float, 2>{grid.points[static_cast<std::size_t>(p)].x, grid.points[static_cast<std::size_t>(p)].y}; };
    for (int k = 0; k < 5; ++k) {
        grid.triangle({0, k + 1, k + 2}, {uv(0), uv(k + 1), uv(k + 2)});
    }
    auto jump = [](const Adjacency& a, int e) {
        const int t = a.edgeCorners(e)[0] / 3;
        return t == 1 || t == 2 ? 1 : -1;
    };
    parts::Part part = tracedPart(grid, {{0, 2}, {0, 4}}, jump, {}, {0});
    extract::Trace trace(part.mesh.faceCount());
    extract::TraceFrame frame(part, trace);
    extract::TraceNode node;
    node.kind = extract::node_kind::point;
    node.reference = 0;
    EXPECT(frame.directionSlot(0, node, 0, 2), std::optional<int>(3));
}

// The length limits of a walk: up to ten nodes always; beyond that given up
// when more than ten steps moved, when the node just added was met before,
// or beyond 200 nodes.
void traceWalkLimits()
{
    extract::FaceLoop loop;
    for (int k = 0; k < 11; ++k) {
        loop.nodes.push_back(k);
    }
    EXPECT(extract::FaceWalker::tooLong(loop, 1), false);
    EXPECT(extract::FaceWalker::tooLong(loop, 0), true);
    loop.nodes.back() = 3;
    EXPECT(extract::FaceWalker::tooLong(loop, 1), true);
    loop.nodes.clear();
    for (int k = 0; k < 200; ++k) {
        loop.nodes.push_back(k);
    }
    EXPECT(extract::FaceWalker::tooLong(loop, 200), false);
    loop.nodes.push_back(200);
    EXPECT(extract::FaceWalker::tooLong(loop, 201), true);
}

// The square of traceSquare with its only start point at (0.75, 0), whose
// coordinates are not whole: the first pass walks nothing, the border points
// with whole coordinates (the four corners) are added and the first pass
// goes with its node. 17 nodes: the 16 grid points, the corners as start
// nodes on their points, and the start node at (0.75, 0) once; the nine
// unit cells are the faces.
void traceReseed()
{
    GridMesh grid;
    for (int j = 0; j < 5; ++j) {
        for (int i = 0; i < 5; ++i) {
            grid.points.push_back({0.75f * static_cast<float>(i), 0.75f * static_cast<float>(j), 0.0f});
        }
    }
    auto uv = [&](int p) { return std::array<float, 2>{grid.points[static_cast<std::size_t>(p)].x, grid.points[static_cast<std::size_t>(p)].y}; };
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            const int a = 5 * j + i;
            grid.triangle({a, a + 1, a + 6}, {uv(a), uv(a + 1), uv(a + 6)});
            grid.triangle({a, a + 6, a + 5}, {uv(a), uv(a + 6), uv(a + 5)});
        }
    }
    parts::Part part = tracedPart(grid, {}, [](const Adjacency&, int) { return 0; }, {}, {1});
    const TracedGrid traced = traceGrid(part, [](Vec3f p) {
        const float x = std::round(p.x);
        const float y = std::round(p.y);
        const bool on = std::fabs(p.x - x) < 1e-5f && std::fabs(p.y - y) < 1e-5f && p.z == 0.0f;
        return on ? static_cast<int>(x) + 4 * static_cast<int>(y) : -1;
    });
    EXPECT(traced.nodes, 17);
    EXPECT(traced.misplaced, 1);
    EXPECT(traced.kinds, (std::array<int, 3>{2, 5, 10}));
    std::vector<std::array<int, 4>> cells;
    for (int y = 0; y < 3; ++y) {
        for (int x = 0; x < 3; ++x) {
            const int k = x + 4 * y;
            cells.push_back({k, k + 1, k + 5, k + 4});
        }
    }
    EXPECT(traced.faces, sortedCells(cells));
}

// The place of a node in a triangle whose corner 1 has no number for (u, v):
// the side 2-0 is the only length that is a number and counts as the
// longest, the node lies on it at the share of its (u, v) distance from
// corner 0.
void tracePlaceWithoutNumber()
{
    const float nan = std::numeric_limits<float>::quiet_NaN();
    GridMesh grid;
    grid.points = {{0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}};
    grid.triangle({0, 1, 2}, {std::array<float, 2>{0.0f, 0.0f}, std::array<float, 2>{nan, nan}, std::array<float, 2>{0.0f, 2.0f}});
    parts::Part part = tracedPart(grid, {}, [](const Adjacency&, int) { return 0; }, {}, {0});
    extract::Trace trace(part.mesh.faceCount());
    const extract::TraceFrame frame(part, trace);
    Vec3f out;
    EXPECT(extract::SegmentTracer::place(frame, 0, -1, {0.0f, 1.0f}, out), true);
    EXPECT(out, (Vec3f{0.0f, 0.5f, 0.0f}));
}

} // namespace

namespace {

std::vector<int> faceStartOf(const Mesh& m)
{
    std::vector<int> starts{0};
    for (int f = 0; f < m.faceCount(); ++f) {
        starts.push_back(m.topology().faceEnd(f));
    }
    return starts;
}

const std::vector<int>& faceGroupValues(const Mesh& m)
{
    return m.maps().findAs<TypedMap<int>>(map_names::faceGroups, ElementKind::Face)->values();
}

// Every rule of the cleanup right after the raw mesh: too few corners, the
// two quad rules, repeats, the first corner against the last one with a
// point, corners without point, a repeat across a corner without point, a
// face left with two corners; unused points stay.
void cleanupFirstPass()
{
    Mesh m = polygonMesh(grid(8), {0, 2, 6, 10, 14, 19, 24, 29, 32, 35},
                         {0, 1,                // too few corners
                          0, 1, 0, 2,          // corners 0 and 2 equal
                          0, 1, 2, 1,          // corners 1 and 3 equal
                          0, 1, 2, 3,          // stays
                          3, 4, 5, 5, 3,       // -> 4 5 3
                          4, 5, 6, 4, -1,      // last point 4 equals the first -> 5 6 4
                          5, -1, 5, 6, 7,      // the repeat across -1 stays -> 5 5 6 7
                          6, 6, 7,             // two corners left
                          2, 3, 4},            // stays
                         FaceMode::Polygons);
    const auto result = DegenerateCleanup::run(m, {});
    EXPECT(result.removedFaces, 4);
    EXPECT(m.topology().mode(), FaceMode::Polygons);
    EXPECT(faceStartOf(m), (std::vector<int>{0, 4, 7, 10, 14, 17}));
    EXPECT(m.topology().cornerVertex(), (std::vector<int>{0, 1, 2, 3, 4, 5, 3, 5, 6, 4, 5, 5, 6, 7, 2, 3, 4}));
    EXPECT(m.pointCount(), 8);

    // Nothing to clean: the quad mesh stays as it is.
    Mesh quads = polygonMesh(grid(6), {0, 4, 8}, {0, 1, 4, 3, 1, 2, 5, 4}, FaceMode::Quads);
    EXPECT(DegenerateCleanup::run(quads, {}).removedFaces, 0);
    EXPECT(quads.topology().mode(), FaceMode::Quads);
    EXPECT(quads.topology().cornerVertex(), (std::vector<int>{0, 1, 4, 3, 1, 2, 5, 4}));
}

// With a face map and given normals: face origin, corner origin of the corner
// actually left out; the unused point removed afterwards.
void cleanupWithOrigins()
{
    Mesh m = polygonMesh(grid(7), {0, 4, 6, 9}, {0, 1, 1, 2, 0, 1, 3, 4, 5}, FaceMode::Polygons);
    m.maps().set(std::make_unique<TypedMap<int>>(std::string(map_names::faceGroups), ElementKind::Face,
                                                 std::vector<int>{10, 11, 12}));
    CustomNormals normals;
    normals.cornerNormal = {100, 101, 102, 103, 104, 105, 106, 107, 108};
    normals.normals.assign(109, Vec3f{0, 0, 1});
    m.customNormals() = normals;
    DegenerateCleanup::Switches switches;
    switches.faceOrigin = true;
    switches.cornerOrigin = true;
    const auto result = DegenerateCleanup::run(m, switches);
    EXPECT(result.removedFaces, 1);
    EXPECT(prep::MeshRepair::removeUnusedPoints(m), 1);
    EXPECT(m.pointCount(), 6);
    EXPECT(faceStartOf(m), (std::vector<int>{0, 3, 6}));
    EXPECT(m.topology().cornerVertex(), (std::vector<int>{0, 1, 2, 3, 4, 5}));
    EXPECT(faceGroupValues(m), (std::vector<int>{10, 12}));
    EXPECT(m.customNormals().has_value(), true);
    if (m.customNormals()) {
        EXPECT(m.customNormals()->cornerNormal, (std::vector<int>{100, 101, 103, 106, 107, 108}));
    }
}

// With the maps kept: a repeated corner goes, the point map stays as it is
// and the hard edge 1-4 follows its half edge.
void cleanupKeepingMaps()
{
    Mesh m = polygonMesh(grid(6), {0, 5, 9}, {0, 1, 1, 4, 3, 1, 2, 5, 4}, FaceMode::Polygons);
    m.maps().set(std::make_unique<TypedMap<float>>(std::string(map_names::density), ElementKind::Point,
                                                   std::vector<float>{0, 1, 2, 3, 4, 5}));
    setSharp(m, {{1, 4}});
    EXPECT(DegenerateCleanup::run(m, {.keepMaps = true}).removedFaces, 0);
    EXPECT(faceStartOf(m), (std::vector<int>{0, 4, 8}));
    EXPECT(m.topology().cornerVertex(), (std::vector<int>{0, 1, 4, 3, 1, 2, 5, 4}));
    const auto* density = m.maps().findAs<TypedMap<float>>(map_names::density, ElementKind::Point);
    EXPECT((density != nullptr && density->values() == std::vector<float>{0, 1, 2, 3, 4, 5}), true);
    EXPECT(isSharp(m, 1, 4), true);
    EXPECT(isSharp(m, 0, 1), false);
}

// Square (0,0) (2,0) (2,2) (0,2) split by a doublet point (1,1).
std::vector<Vec3f> doubletPoints()
{
    return {{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {1, 1, 0}};
}

// The doublet at the border: both faces become one polygon from the first
// end point of the smallest border edge, in the direction of the first face;
// the doublet point stays unused, the corners at 90 degrees stay.
void doubletAtBorder()
{
    Mesh m = polygonMesh(doubletPoints(), {0, 4, 8}, {0, 1, 2, 4, 0, 4, 2, 3}, FaceMode::Quads);
    EXPECT(extract::Doublets::run(m), 1);
    EXPECT(m.topology().mode(), FaceMode::Polygons);
    EXPECT(faceStartOf(m), (std::vector<int>{0, 4}));
    EXPECT(m.topology().cornerVertex(), (std::vector<int>{0, 1, 2, 3}));
    EXPECT(m.pointCount(), 5);
    EXPECT(m.adjacency().pointCorners(4).empty(), true);
}

// The same doublet inside a frame of four quads: the frame keeps its order,
// the new quad comes last.
void doubletInside()
{
    std::vector<Vec3f> p = doubletPoints();
    p.insert(p.end(), {{-1, -1, 0}, {3, -1, 0}, {3, 3, 0}, {-1, 3, 0}});
    Mesh m = polygonMesh(p, {0, 4, 8, 12, 16, 20, 24},
                         {0, 1, 2, 4, 0, 4, 2, 3, 5, 6, 1, 0, 6, 7, 2, 1, 7, 8, 3, 2, 8, 5, 0, 3}, FaceMode::Quads);
    EXPECT(extract::Doublets::run(m), 1);
    EXPECT(faceStartOf(m), (std::vector<int>{0, 4, 8, 12, 16, 20}));
    EXPECT(m.topology().cornerVertex(),
           (std::vector<int>{5, 6, 1, 0, 6, 7, 2, 1, 7, 8, 3, 2, 8, 5, 0, 3, 0, 1, 2, 3}));
    EXPECT(m.adjacency().nonManifoldEdges().empty(), true);
}

// A doublet (1,1) whose dissolving leaves the point (1,0) between two
// straight edges: it is left out, and the face below it is built anew as a
// group of its own (new start corner, appended); the face map follows the
// smallest face of each group. Without straight points the face below stays.
void dissolveStraightPoint()
{
    const std::vector<Vec3f> p = {{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {1, 1, 0}, {1, 0, 0}, {1, -1, 0}};
    const std::vector<int> starts = {0, 4, 9, 13};
    const std::vector<int> corners = {5, 1, 2, 4, 5, 4, 2, 3, 0, 5, 0, 6, 1};
    auto withGroups = [&]() {
        Mesh m = polygonMesh(p, starts, corners, FaceMode::Polygons);
        m.maps().set(std::make_unique<TypedMap<int>>(std::string(map_names::faceGroups), ElementKind::Face,
                                                     std::vector<int>{10, 11, 12}));
        return m;
    };

    Mesh m = withGroups();
    EXPECT(extract::Doublets::run(m), 1);
    EXPECT(faceStartOf(m), (std::vector<int>{0, 4, 7}));
    EXPECT(m.topology().cornerVertex(), (std::vector<int>{0, 1, 2, 3, 0, 6, 1}));
    EXPECT(faceGroupValues(m), (std::vector<int>{10, 12}));
    EXPECT(m.adjacency().nonManifoldEdges().empty(), true);

    Mesh kept = withGroups();
    const auto doubletEdges = kept.adjacency().pointEdges(4);
    const std::vector<int> edges(doubletEdges.begin(), doubletEdges.end());
    extract::EdgeDissolve::run(kept, edges, false);
    EXPECT(faceStartOf(kept), (std::vector<int>{0, 4, 9}));
    EXPECT(kept.topology().cornerVertex(), (std::vector<int>{5, 0, 6, 1, 0, 5, 1, 2, 3}));
    EXPECT(faceGroupValues(kept), (std::vector<int>{12, 10}));
}

// A ring of four quads dissolved into one group: every border loop becomes a
// polygon of its own, the inner one against the outer direction.
void dissolveRing()
{
    const std::vector<Vec3f> p = {{0, 0, 0}, {3, 0, 0}, {3, 3, 0}, {0, 3, 0},
                                  {1, 1, 0}, {2, 1, 0}, {2, 2, 0}, {1, 2, 0}};
    Mesh m = polygonMesh(p, {0, 4, 8, 12, 16}, {0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7}, FaceMode::Quads);
    const Adjacency& a = m.adjacency();
    const std::vector<int> edges = {a.edgeBetween(0, 4), a.edgeBetween(1, 5), a.edgeBetween(2, 6), a.edgeBetween(3, 7)};
    extract::EdgeDissolve::run(m, edges, true);
    EXPECT(faceStartOf(m), (std::vector<int>{0, 4, 8}));
    EXPECT(m.topology().cornerVertex(), (std::vector<int>{0, 1, 2, 3, 4, 7, 6, 5}));
}

extract::OutputPoints outputPoints(int count)
{
    extract::OutputPoints points;
    points.onBorderEdge.resize(static_cast<std::size_t>(count));
    points.reachedBorder.resize(static_cast<std::size_t>(count));
    points.source.resize(static_cast<std::size_t>(count));
    points.weight.assign(static_cast<std::size_t>(count), 1.0f);
    return points;
}

// Triangle of area 1000.
Mesh areaTriangle()
{
    return polygonMesh({{0, 0, 0}, {100, 0, 0}, {0, 20, 0}}, {0, 3}, {0, 1, 2}, FaceMode::Triangles);
}

// Three triangles on one edge: both of its end points get a point per
// further fan; the point data follow the source point.
void holeSplitFans()
{
    Mesh m = polygonMesh(grid(5), {0, 3, 6, 9}, {0, 1, 2, 1, 0, 3, 0, 1, 4}, FaceMode::Triangles);
    extract::OutputPoints points = outputPoints(5);
    points.onBorderEdge.set(0);
    points.reachedBorder.set(1);
    points.source = {{1, 10}, {2, 11}, {4, 12}, {0, -1}, {1, 14}};
    points.weight = {0.5f, 0.25f, 1.0f, 1.0f, 0.0f};
    points.triangle = {0, 1, 2, 3, 4};
    for (auto& n : points.neighbours) {
        n.assign(5, 9);
    }
    const Mesh original = areaTriangle();
    extract::HoleClosing(m, points, original, 1, nullptr, {.rebuildPolygons = false, .followBorders = false}).run();
    EXPECT(m.pointCount(), 9);
    EXPECT(m.topology().cornerVertex(), (std::vector<int>{0, 1, 2, 7, 5, 3, 6, 8, 4}));
    EXPECT(m.points()[5] == m.points()[0] && m.points()[8] == m.points()[1], true);
    std::vector<int> onBorder;
    std::vector<int> reached;
    for (std::size_t v = 0; v < 9; ++v) {
        if (points.onBorderEdge.testBounded(v)) {
            onBorder.push_back(static_cast<int>(v));
        }
        if (points.reachedBorder.testBounded(v)) {
            reached.push_back(static_cast<int>(v));
        }
    }
    EXPECT(onBorder, (std::vector<int>{0, 5, 6}));
    EXPECT(reached, (std::vector<int>{1, 7, 8}));
    EXPECT(points.source, (std::vector<extract::PointSource>{{1, 10}, {2, 11}, {4, 12}, {0, -1}, {1, 14},
                                                             {1, 10}, {1, 10}, {2, 11}, {2, 11}}));
    EXPECT(points.weight, (std::vector<float>{0.5f, 0.25f, 1.0f, 1.0f, 0.0f, 0.5f, 0.5f, 0.25f, 0.25f}));
    EXPECT(points.triangle, (std::vector<int>{0, 1, 2, 3, 4, -1, -1, -1, -1}));
    EXPECT(points.neighbours[3], (std::vector<int>{9, 9, 9, 9, 9, -1, -1, -1, -1}));
}

// Ring of four quads around a hole; the outer points reached the border.
Mesh holeRing(const std::vector<Vec3f>& inner)
{
    std::vector<Vec3f> p = {{0, 0, 0}, {3, 0, 0}, {3, 3, 0}, {0, 3, 0}};
    p.insert(p.end(), inner.begin(), inner.end());
    return polygonMesh(p, {0, 4, 8, 12, 16}, {0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7}, FaceMode::Polygons);
}

extract::OutputPoints ringPoints()
{
    extract::OutputPoints points = outputPoints(8);
    for (std::size_t v = 0; v < 4; ++v) {
        points.reachedBorder.set(v);
    }
    return points;
}

// The inner loop (5 4 7 6) closes against its direction; a thin hole closes
// by joining its short diagonal (q0 to q2, turned by a quarter q1 to q3)
// while it is below area / (500 budget); the outer loop reached the border
// and stays open.
void holesAndCollapse()
{
    const Mesh original = areaTriangle();
    const extract::HoleClosing::Options noNotches{.rebuildPolygons = true, .followBorders = false};
    const std::vector<int> ringCorners = {0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7};

    Mesh square = holeRing({{1, 1, 0}, {2, 1, 0}, {2, 2, 0}, {1, 2, 0}});
    extract::OutputPoints points = ringPoints();
    extract::HoleClosing(square, points, original, 1, nullptr, noNotches).run();
    EXPECT(faceStartOf(square), (std::vector<int>{0, 4, 8, 12, 16, 20}));
    std::vector<int> expected = ringCorners;
    expected.insert(expected.end(), {6, 7, 4, 5});
    EXPECT(square.topology().cornerVertex(), expected);

    const std::vector<Vec3f> thin = {{1.5f, 1.4f, 0}, {2, 1.5f, 0}, {1.5f, 1.6f, 0}, {1, 1.5f, 0}};
    Mesh joined = holeRing(thin);
    points = ringPoints();
    extract::HoleClosing(joined, points, original, 1, nullptr, noNotches).run();
    EXPECT(faceStartOf(joined), (std::vector<int>{0, 4, 8, 12, 16}));
    EXPECT(joined.topology().cornerVertex(), (std::vector<int>{0, 1, 5, 4, 1, 2, 4, 5, 2, 3, 7, 4, 3, 0, 4, 7}));

    // The same hole turned by a quarter: the short diagonal joins q1 to q3.
    const std::vector<Vec3f> turned = {{1.5f, 1, 0}, {1.6f, 1.5f, 0}, {1.5f, 2, 0}, {1.4f, 1.5f, 0}};
    Mesh across = holeRing(turned);
    points = ringPoints();
    extract::HoleClosing(across, points, original, 1, nullptr, noNotches).run();
    EXPECT(faceStartOf(across), (std::vector<int>{0, 4, 8, 12, 16}));
    EXPECT(across.topology().cornerVertex(), (std::vector<int>{0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 5, 6, 3, 0, 4, 5}));

    Mesh large = holeRing(thin);
    points = ringPoints();
    extract::HoleClosing(large, points, original, 1000, nullptr, noNotches).run();
    EXPECT(large.topology().cornerVertex(), expected);
}

// Without frozen borders: the outer loop is cut at its points on border
// edges (0 and 2) into polygons against the loop direction, after the hole.
void holeNotches()
{
    const Mesh original = areaTriangle();
    Mesh m = holeRing({{1, 1, 0}, {2, 1, 0}, {2, 2, 0}, {1, 2, 0}});
    extract::OutputPoints points = ringPoints();
    points.onBorderEdge.set(0);
    points.onBorderEdge.set(2);
    extract::HoleClosing(m, points, original, 1, nullptr, {}).run();
    EXPECT(faceStartOf(m), (std::vector<int>{0, 4, 8, 12, 16, 20, 23, 26}));
    EXPECT(std::vector<int>(m.topology().cornerVertex().begin() + 16, m.topology().cornerVertex().end()),
           (std::vector<int>{6, 7, 4, 5, 2, 1, 0, 0, 3, 2}));
}

// Frozen borders: the ring between the border loop (0 4 5 6 1 2 3) and the
// kept chain (3 2 1 0) is cut at the marked points, round to the last mark
// again; one mark gives the whole ring from it.
void holeGapRings()
{
    const Mesh original = polygonMesh({{0, 0, 0}, {4, 0, 0}, {4, 4, 0}, {0, 4, 0}}, {0, 3, 6}, {0, 1, 2, 0, 2, 3},
                                      FaceMode::Triangles);
    const std::vector<FreezeMark> marks(4, FreezeMark{1});
    const std::vector<int> originalToOutput = {0, 1, 2, 3};
    const std::vector<Vec3f> p = {{0, 0, 0}, {4, 0, 0}, {4, 4, 0}, {0, 4, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}};
    auto run = [&](const BitArray& locked, extract::PointSource middle) {
        Mesh m = polygonMesh(p, {0, 7}, {0, 4, 5, 6, 1, 2, 3}, FaceMode::Polygons);
        extract::OutputPoints points = outputPoints(7);
        for (std::size_t v = 0; v < 7; ++v) {
            points.reachedBorder.set(v);
        }
        points.source = {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {4, 0}, middle, {4, 0}};
        const extract::FrozenBorderInput frozen{marks, originalToOutput, &locked};
        extract::HoleClosing(m, points, original, 1, &frozen, {}).run();
        return std::pair{faceStartOf(m), m.topology().cornerVertex()};
    };
    const int bottom = original.adjacency().edgeBetween(0, 1);
    BitArray lockBottom(static_cast<std::size_t>(original.adjacency().edgeCount()));
    lockBottom.set(static_cast<std::size_t>(bottom));
    const auto [cutStarts, cutCorners] = run(lockBottom, {2, bottom});
    EXPECT(cutStarts, (std::vector<int>{0, 7, 10, 13}));
    EXPECT(cutCorners, (std::vector<int>{0, 4, 5, 6, 1, 2, 3, 1, 6, 5, 5, 4, 0}));

    BitArray lockRight(static_cast<std::size_t>(original.adjacency().edgeCount()));
    lockRight.set(static_cast<std::size_t>(original.adjacency().edgeBetween(1, 2)));
    const auto [oneStarts, oneCorners] = run(lockRight, {2, bottom});
    EXPECT(oneStarts, (std::vector<int>{0, 7, 12}));
    EXPECT(oneCorners, (std::vector<int>{0, 4, 5, 6, 1, 2, 3, 1, 6, 5, 4, 0}));
}

} // namespace

namespace {

// Number of faces by corner count (index 7: seven and more).
std::vector<int> cornerCounts(const Mesh& m)
{
    std::vector<int> counts(8, 0);
    for (int f = 0; f < m.faceCount(); ++f) {
        ++counts[static_cast<std::size_t>(std::min(m.topology().faceSize(f), 7))];
    }
    return counts;
}

// Number of output points by source: none, point, edge, triangle.
std::array<int, 4> sourceKinds(const extract::OutputPoints& points)
{
    std::array<int, 4> out{};
    for (const extract::PointSource& s : points.source) {
        ++out[static_cast<std::size_t>(s.type == 4 ? 3 : s.type)];
    }
    return out;
}

// A loop of the given nodes for the face closing, every side marked as in `degenerate`.
extract::FaceLoop faceLoop(std::vector<int> nodes, std::vector<std::uint8_t> degenerate)
{
    extract::FaceLoop loop;
    loop.directions.assign(nodes.size(), 0);
    loop.triangles.assign(nodes.size(), 0);
    loop.nodes = std::move(nodes);
    loop.degenerate = std::move(degenerate);
    return loop;
}

// A trace with `count` free nodes in one triangle, numbered in order.
extract::Trace freeNodes(int count, int kind = extract::node_kind::free)
{
    extract::Trace trace(1);
    for (int i = 0; i < count; ++i) {
        const int k = trace.nodes.create(0, {static_cast<float>(i), 0.0f, 0.0f}, {0.0f, 0.0f}, false);
        trace.nodes[k].id = trace.nodes.takeId();
        trace.nodes[k].kind = kind;
        trace.nodes[k].reference = i;
    }
    return trace;
}

std::vector<std::vector<int>> polygonRows(const Csr& polygons)
{
    std::vector<std::vector<int>> rows;
    for (int r = 0; r < polygons.rows(); ++r) {
        rows.push_back(v(polygons.row(r)));
    }
    return rows;
}

// Closing faces: a spike in a heptagon gives two quads and a triangle; a
// degenerate side keeps the node with the smaller number; a triangle without
// node on a mesh point that keeps a degenerate side (two crease chains) is
// dropped.
void faceClosingRules()
{
    extract::FaceClosing closer(nullptr);
    {
        extract::Trace trace = freeNodes(6);
        extract::FaceLoop loop = faceLoop({0, 1, 2, 1, 3, 4, 5}, std::vector<std::uint8_t>(7, 0));
        closer.close(loop, trace);
        EXPECT(polygonRows(trace.polygons), (std::vector<std::vector<int>>{{2, 1, 3, 4}, {5, 0, 1, 2}, {2, 4, 5}}));
    }
    {
        extract::Trace trace = freeNodes(5);
        extract::FaceLoop loop = faceLoop({0, 1, 4, 3, 2}, {0, 0, 1, 0, 0});
        closer.close(loop, trace);
        EXPECT(polygonRows(trace.polygons), (std::vector<std::vector<int>>{{0, 1, 3, 2}}));
        EXPECT(trace.nodes[4].mergeTarget, 3);
    }
    {
        extract::Trace trace = freeNodes(3, extract::node_kind::edge);
        trace.curves.edges = BitArray(3, true);
        trace.curves.edgeChain = {0, 1, 1};
        extract::FaceLoop loop = faceLoop({0, 1, 2}, {1, 0, 0});
        closer.close(loop, trace);
        EXPECT(trace.polygons.rows(), 0);
    }
}

// Special ways of closing faces. A walk from a seed that comes back over its
// second node keeps that one (the seed has the smaller number). With frozen
// borders, nodes on points of the triangle mesh: two kept points stay apart;
// a kept point stays against a point that is not kept (also with the larger
// number), and so does the other one when it is the kept one; a kept point
// and a node on a crease edge stay apart, both ways round. A hexagon with a
// spike loses the spike and the node met again. A small face that doubles
// back on itself only by its numbers (a merged node) is left out.
void faceClosingBranches()
{
    auto closed = [](const std::vector<FreezeMark>* marks, extract::Trace& trace, extract::FaceLoop loop) {
        extract::FaceClosing closer(marks);
        closer.close(loop, trace);
        return polygonRows(trace.polygons);
    };
    using Rows = std::vector<std::vector<int>>;
    {
        extract::Trace trace = freeNodes(4);
        trace.nodes[0].flags |= extract::node_flag::seed;
        EXPECT(closed(nullptr, trace, faceLoop({0, 1, 2, 3, 1}, {1, 0, 0, 0, 0})), (Rows{{1, 2, 3, 1}}));
        EXPECT(trace.nodes[0].mergeTarget, 1);
        EXPECT(trace.nodes[1].mergeTarget, -1);
    }
    auto marks = [](std::initializer_list<int> kept) {
        std::vector<FreezeMark> out(4);
        for (const int p : kept) {
            out[static_cast<std::size_t>(p)].value = 1;
        }
        return out;
    };
    {
        extract::Trace trace = freeNodes(4, extract::node_kind::point);
        const std::vector<FreezeMark> both = marks({0, 1});
        EXPECT(closed(&both, trace, faceLoop({0, 1, 2, 3}, {1, 0, 0, 0})), (Rows{{0, 1, 2, 3}}));
        EXPECT(trace.nodes[0].mergeTarget, -1);
        EXPECT(trace.nodes[1].mergeTarget, -1);
    }
    {
        extract::Trace trace = freeNodes(4, extract::node_kind::point);
        const std::vector<FreezeMark> first = marks({1});
        EXPECT(closed(&first, trace, faceLoop({1, 0, 2, 3}, {1, 0, 0, 0})), (Rows{{1, 2, 3}}));
        EXPECT(trace.nodes[0].mergeTarget, 1);
    }
    {
        extract::Trace trace = freeNodes(4, extract::node_kind::point);
        const std::vector<FreezeMark> second = marks({1});
        EXPECT(closed(&second, trace, faceLoop({0, 1, 2, 3}, {1, 0, 0, 0})), (Rows{{1, 2, 3}}));
        EXPECT(trace.nodes[0].mergeTarget, 1);
    }
    for (const int onCrease : {1, 0}) {
        extract::Trace trace = freeNodes(4, extract::node_kind::point);
        trace.nodes[onCrease].kind = extract::node_kind::edge;
        trace.curves.edges = BitArray(4);
        trace.curves.edges.set(static_cast<std::size_t>(onCrease));
        const std::vector<FreezeMark> kept = marks({1 - onCrease});
        EXPECT(closed(&kept, trace, faceLoop({0, 1, 2, 3}, {1, 0, 0, 0})), (Rows{{0, 1, 2, 3}}));
        EXPECT(trace.nodes[0].mergeTarget, -1);
        EXPECT(trace.nodes[1].mergeTarget, -1);
    }
    {
        extract::Trace trace = freeNodes(5);
        EXPECT(closed(nullptr, trace, faceLoop({0, 1, 2, 1, 3, 4}, std::vector<std::uint8_t>(6, 0))),
               (Rows{{0, 1, 3, 4}}));
    }
    {
        extract::Trace trace = freeNodes(4);
        trace.nodes[3].mergeTarget = 1;
        EXPECT(closed(nullptr, trace, faceLoop({0, 1, 2, 3}, std::vector<std::uint8_t>(4, 0))).size(),
               std::size_t{0});
    }
}

// Cutting an octagon around a 3 x 1 rectangle into quads: the shortest
// diagonal first (corners 2 and 5), then, with the positions cut as well,
// corners 1 and 6; the rest 0 1 6 7. Without positions nothing is cut.
void splitIntoQuads()
{
    std::vector<int> numbers{0, 1, 2, 3, 4, 5, 6, 7};
    const std::vector<Vec3f> positions{{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0},
                                       {3, 1, 0}, {2, 1, 0}, {1, 1, 0}, {0, 1, 0}};
    const auto quads = extract::FaceClosing::splitIntoQuads(numbers, positions);
    EXPECT(quads, (std::vector<std::array<int, 4>>{{2, 3, 4, 5}, {1, 2, 5, 6}}));
    EXPECT(numbers, (std::vector<int>{0, 1, 6, 7}));
    const float nan = std::numeric_limits<float>::quiet_NaN();
    std::vector<int> undecided{0, 1, 2, 3, 4};
    EXPECT(extract::FaceClosing::splitIntoQuads(undecided, std::vector<Vec3f>(5, Vec3f{nan, nan, nan})).size(),
           std::size_t{0});
    EXPECT(undecided.size(), std::size_t{5});
}

// Border classification on a raw heptagon over a strip of two squares: the
// point between the anchors on points 0 and 1 lies 0.00005 off the border
// edge and is marked, the one between points 1 and 2 lies 0.001 off and is
// not. The loop is looked at from its first anchor round and back to it.
void borderClassification()
{
    parts::Part part;
    part.mesh = polygonMesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}}, {0, 3, 6, 9, 12},
                            {0, 1, 4, 0, 4, 3, 1, 2, 5, 1, 5, 4}, FaceMode::Triangles);
    extract::RawMesh raw;
    raw.mesh = polygonMesh({{0, 0, 0}, {0.5f, 0.00005f, 0}, {1, 0, 0}, {1.5f, 0.001f, 0}, {2, 0, 0}, {2, 1, 0}, {0, 1, 0}},
                           {0, 7}, {0, 1, 2, 3, 4, 5, 6}, FaceMode::Polygons);
    raw.points.source = {{1, 0}, {4, 0}, {1, 1}, {4, 2}, {1, 2}, {1, 5}, {1, 3}};
    raw.points.onBorderEdge = BitArray(7);
    for (const int p : {0, 2, 4, 5, 6}) {
        raw.points.onBorderEdge.set(static_cast<std::size_t>(p));
    }
    EXPECT(extract::BorderClassifier::visitOrder(std::vector<int>{1, 0, 3, 2}, raw), (std::vector<int>{0, 3, 2, 1, 0}));
    extract::BorderClassifier::run(raw, part);
    EXPECT(raw.points.onBorderEdge.test(1), true);
    EXPECT(raw.points.onBorderEdge.test(3), false);
    EXPECT(raw.mesh.faceCount(), 1);
}

// The square [0, 3]^2 of the trace check: lattice of step 0.75, cells split
// along the rising diagonal, (u, v) = (x, y).
GridMesh squareLattice()
{
    GridMesh grid;
    for (int j = 0; j < 5; ++j) {
        for (int i = 0; i < 5; ++i) {
            grid.points.push_back({0.75f * static_cast<float>(i), 0.75f * static_cast<float>(j), 0.0f});
        }
    }
    auto uv = [&](int p) {
        return std::array<float, 2>{grid.points[static_cast<std::size_t>(p)].x, grid.points[static_cast<std::size_t>(p)].y};
    };
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            const int a = 5 * j + i;
            grid.triangle({a, a + 1, a + 6}, {uv(a), uv(a + 1), uv(a + 6)});
            grid.triangle({a, a + 6, a + 5}, {uv(a), uv(a + 6), uv(a + 5)});
        }
    }
    return grid;
}

void rawMeshOf(parts::Part& part)
{
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    extract::Extraction::rawMesh(part, {}, progress);
}

// Raw mesh of the square: the nine cells as quads, manifold, one border
// loop; the three corners besides the start point lie on mesh points and
// are snapped onto them (source point, place of the point), the eight
// border nodes and the two on the diagonals come from edges, the two others
// from triangles. Only the start node was on a point before the snapping.
void rawSquare()
{
    parts::Part part = tracedPart(squareLattice(), {}, [](const Adjacency&, int) { return 0; }, {}, {0});
    part.quadBudget = 9;
    rawMeshOf(part);
    EXPECT(part.raw.has_value(), true);
    if (!part.raw) {
        return;
    }
    const extract::RawMesh& raw = *part.raw;
    EXPECT(cornerCounts(raw.mesh), (std::vector<int>{0, 0, 0, 0, 9, 0, 0, 0}));
    EXPECT(raw.mesh.adjacency().nonManifoldEdges().size(), std::size_t{0});
    EXPECT(BoundaryLoops(raw.mesh).count(), 1);
    EXPECT(raw.mesh.pointCount(), 16);
    EXPECT(sourceKinds(raw.points), (std::array<int, 4>{0, 4, 10, 2}));
    std::vector<int> onPoints;
    int placed = 0;
    for (std::size_t p = 0; p < raw.points.source.size(); ++p) {
        const extract::PointSource s = raw.points.source[p];
        if (s.type == 1) {
            onPoints.push_back(s.id);
            placed += raw.mesh.points()[p] == part.mesh.points()[static_cast<std::size_t>(s.id)] ? 1 : 0;
        }
    }
    std::sort(onPoints.begin(), onPoints.end());
    EXPECT(onPoints, (std::vector<int>{0, 4, 20, 24}));
    EXPECT(placed, 4);
    EXPECT(std::count_if(raw.originalToOutput.begin(), raw.originalToOutput.end(), [](int o) { return o >= 0; }), 1);
    EXPECT(raw.points.onBorderEdge.count(), std::size_t{12});
}

// Raw mesh of the plate around a point with three quads (the lattice of the
// trace check): twelve quads, manifold, one border loop of twelve points,
// the centre the only inner point without four edges; the spoke ends and the
// outer corners of the sectors lie on mesh points and are snapped.
void rawStar()
{
    const std::array<Vec3f, 3> spoke{Vec3f{0.0f, 1.0f, 0.0f}, Vec3f{-0.8660254f, -0.5f, 0.0f},
                                     Vec3f{0.8660254f, -0.5f, 0.0f}};
    auto at = [&](int i, float s, float t) {
        const Vec3f a = spoke[static_cast<std::size_t>(i)];
        const Vec3f b = spoke[static_cast<std::size_t>((i + 1) % 3)];
        return Vec3f{s * a.x + t * b.x, s * a.y + t * b.y, s * a.z + t * b.z};
    };
    auto point = [](int i, int si, int ti) {
        if (si == 0 && ti == 0) {
            return 0;
        }
        if (ti == 0) {
            return 1 + 3 * i + si - 1;
        }
        if (si == 0) {
            return 1 + 3 * ((i + 1) % 3) + ti - 1;
        }
        return 10 + 9 * i + 3 * (si - 1) + ti - 1;
    };
    auto step = [](int k) { return static_cast<float>(2 * k) / 3.0f; };
    GridMesh grid;
    grid.points.resize(37);
    for (int i = 0; i < 3; ++i) {
        for (int si = 0; si < 4; ++si) {
            for (int ti = 0; ti < 4; ++ti) {
                grid.points[static_cast<std::size_t>(point(i, si, ti))] = at(i, step(si), step(ti));
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        auto uv = [&](int si, int ti) {
            const float s = step(si);
            const float t = step(ti);
            return i == 0 ? std::array<float, 2>{s, t}
                          : (i == 1 ? std::array<float, 2>{0.0f - t, s} : std::array<float, 2>{0.0f - s, 0.0f - t});
        };
        for (int si = 0; si < 3; ++si) {
            for (int ti = 0; ti < 3; ++ti) {
                grid.triangle({point(i, si, ti), point(i, si + 1, ti), point(i, si + 1, ti + 1)},
                              {uv(si, ti), uv(si + 1, ti), uv(si + 1, ti + 1)});
                grid.triangle({point(i, si, ti), point(i, si + 1, ti + 1), point(i, si, ti + 1)},
                              {uv(si, ti), uv(si + 1, ti + 1), uv(si, ti + 1)});
            }
        }
    }
    auto jump = [](const Adjacency& a, int e) { return a.edgeCorners(e)[0] / 3 < 18 ? 1 : -1; };
    parts::Part part = tracedPart(grid, {{0, 1}, {1, 2}, {2, 3}}, jump, {0}, {});
    part.quadBudget = 12;
    rawMeshOf(part);
    EXPECT(part.raw.has_value(), true);
    if (!part.raw) {
        return;
    }
    const extract::RawMesh& raw = *part.raw;
    EXPECT(cornerCounts(raw.mesh), (std::vector<int>{0, 0, 0, 0, 12, 0, 0, 0}));
    EXPECT(raw.mesh.adjacency().nonManifoldEdges().size(), std::size_t{0});
    EXPECT(raw.mesh.adjacency().counts().irregularInnerPoints, 1);
    const BoundaryLoops loops(raw.mesh);
    EXPECT(loops.count(), 1);
    EXPECT(loops.count() == 1 ? loops.points().rowSize(0) : 0, 12);
    EXPECT(sourceKinds(raw.points), (std::array<int, 4>{0, 7, 12, 0}));
}

// Raw mesh of the open cylinder of the trace check: eight quads, two border
// loops of four points, every inner point with four edges; the grid points
// on whole columns of the lattice (u = 0 and 2 on both borders) are snapped.
void rawCylinder()
{
    const float radius = 4.0f / (2.0f * std::numbers::pi_v<float>);
    auto at = [&](float u, float v) {
        const float angle = 2.0f * std::numbers::pi_v<float> * u / 4.0f;
        return Vec3f{radius * std::cos(angle), radius * std::sin(angle), v};
    };
    auto step = [](int k) { return static_cast<float>(2 * k) / 3.0f; };
    GridMesh grid;
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 6; ++col) {
            grid.points.push_back(at(step(col), step(row)));
        }
    }
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 6; ++col) {
            const int next = (col + 1) % 6;
            const std::array<float, 2> a{step(col), step(row)};
            const std::array<float, 2> b{step(col + 1), step(row)};
            const std::array<float, 2> c{step(col + 1), step(row + 1)};
            const std::array<float, 2> d{step(col), step(row + 1)};
            const int pa = 6 * row + col;
            const int pb = 6 * row + next;
            const int pc = 6 * (row + 1) + next;
            const int pd = 6 * (row + 1) + col;
            grid.triangle({pa, pb, pc}, {a, b, c});
            grid.triangle({pa, pc, pd}, {a, c, d});
        }
    }
    parts::Part part = tracedPart(grid, {{0, 6}, {6, 12}, {12, 18}}, [](const Adjacency&, int) { return 0; }, {}, {0});
    part.quadBudget = 8;
    rawMeshOf(part);
    EXPECT(part.raw.has_value(), true);
    if (!part.raw) {
        return;
    }
    const extract::RawMesh& raw = *part.raw;
    EXPECT(cornerCounts(raw.mesh), (std::vector<int>{0, 0, 0, 0, 8, 0, 0, 0}));
    EXPECT(raw.mesh.adjacency().nonManifoldEdges().size(), std::size_t{0});
    EXPECT(raw.mesh.adjacency().counts().irregularInnerPoints, 0);
    const BoundaryLoops loops(raw.mesh);
    EXPECT(loops.count(), 2);
    std::vector<int> loopSizes;
    for (int r = 0; r < loops.count(); ++r) {
        loopSizes.push_back(loops.points().rowSize(r));
    }
    EXPECT(loopSizes, (std::vector<int>{4, 4}));
    EXPECT(sourceKinds(raw.points), (std::array<int, 4>{0, 4, 8, 0}));
}

// The square with frozen borders (all sixteen border points kept): the
// trace takes the kept points between two grid nodes into the faces, so the
// four corner cells become hexagons and the four other border cells
// pentagons, the middle cell stays a quad. The three corners that end on
// mesh points only by snapping get extra output points (the start corner is
// a node on its point). The border classification then dissolves the inner
// edges at the eight border nodes between the kept points: none of them is
// left in a face (the two nodes on the diagonals inside stay).
void rawFrozenPlate()
{
    parts::Part part = tracedPart(squareLattice(), {}, [](const Adjacency&, int) { return 0; }, {}, {0});
    freezeBorders(part, {0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24});
    part.quadBudget = 9;
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    extract::FaceClosing closer(part.frozenBorderMarks());
    extract::IsolineTracer::trace(part, closer, progress);
    EXPECT(part.trace.has_value(), true);
    if (!part.trace) {
        return;
    }
    extract::RawMesh raw = extract::QuadMeshBuilder::build(part, *part.trace);
    EXPECT(cornerCounts(raw.mesh), (std::vector<int>{0, 0, 0, 0, 1, 4, 4, 0}));
    const int numbered = part.trace->nodes.idCount();
    EXPECT(raw.mesh.pointCount() - numbered, 3);
    std::vector<int> extra;
    for (std::size_t p = static_cast<std::size_t>(numbered); p < raw.points.source.size(); ++p) {
        extra.push_back(raw.points.source[p].type == 1 ? raw.points.source[p].id : -1);
    }
    std::sort(extra.begin(), extra.end());
    EXPECT(extra, (std::vector<int>{4, 20, 24}));
    EXPECT(std::count_if(raw.originalToOutput.begin(), raw.originalToOutput.end(), [](int o) { return o >= 0; }), 16);
    extract::BorderClassifier::run(raw, part);
    // (1, 1) and (2, 2) on the diagonals come from edges as well but lie inside.
    int onEdges = 0;
    for (const int p : raw.mesh.topology().cornerVertex()) {
        const auto bit = static_cast<std::size_t>(p);
        onEdges += raw.points.source[bit].type == 2 && raw.points.onBorderEdge.testBounded(bit) ? 1 : 0;
    }
    EXPECT(onEdges, 0);
}

// The square with frozen borders through the gap rings: the three corners
// snapped onto their points lie exactly on their extra points, and both are
// corners of faces after the rings are closed, but never of the same face:
// the snapped node only on the edges of the raw mesh, the extra point only
// on the edges of the kept border chain. No edge has length zero.
void rawFrozenPlateRings()
{
    parts::Part part = tracedPart(squareLattice(), {}, [](const Adjacency&, int) { return 0; }, {}, {0});
    freezeBorders(part, {0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24});
    part.quadBudget = 9;
    rawMeshOf(part);
    EXPECT(part.raw.has_value(), true);
    if (!part.raw) {
        return;
    }
    const Mesh& m = part.raw->mesh;
    const Topology& topology = m.topology();
    auto place = [&](int c) { return m.points()[static_cast<std::size_t>(topology.vertex(c))]; };
    int zeroEdges = 0;
    int samePlaceInFace = 0;
    int samePlaceElsewhere = 0;
    for (int f = 0; f < topology.faceCount(); ++f) {
        for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
            const int next = c + 1 < topology.faceEnd(f) ? c + 1 : topology.faceBegin(f);
            zeroEdges += place(c) == place(next) ? 1 : 0;
            for (int d = c + 1; d < topology.faceEnd(f); ++d) {
                samePlaceInFace += place(c) == place(d) ? 1 : 0;
            }
        }
    }
    const std::vector<int> used = [&] {
        std::vector<int> out = topology.cornerVertex();
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }();
    for (std::size_t i = 0; i < used.size(); ++i) {
        for (std::size_t j = i + 1; j < used.size(); ++j) {
            samePlaceElsewhere += m.points()[static_cast<std::size_t>(used[i])] ==
                                          m.points()[static_cast<std::size_t>(used[j])]
                                      ? 1
                                      : 0;
        }
    }
    EXPECT(zeroEdges, 0);
    EXPECT(samePlaceInFace, 0);
    EXPECT(samePlaceElsewhere, 3);
}

} // namespace

namespace {

// The last cleanup: the quad with corners 0 and 2 on the same point goes,
// the repeated corner of the third face goes, the four points without corner
// go. The colours stay per point that stays, the hard edge 2-3 stays, the
// point data follow the new numbering (a shorter list stays as short over
// the points it covered), and the curve edges are dropped.
void finalCleanupMapsAndUnused()
{
    extract::RawMesh raw;
    raw.mesh = polygonMesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0, 0}, {3, 0, 0}, {2, 1, 0}, {5, 5, 5},
                            {1, 2, 0}},
                           {0, 4, 8, 12}, {0, 1, 2, 3, 4, 5, 4, 6, 3, 2, 8, 8}, FaceMode::Polygons);
    std::vector<Color4ub> colors;
    for (int i = 0; i < 9; ++i) {
        colors.push_back({static_cast<std::uint8_t>(10 * i), 1, 2, 255});
    }
    raw.mesh.maps().set(std::make_unique<TypedMap<Color4ub>>(std::string(map_names::colors), ElementKind::Point, colors));
    setSharp(raw.mesh, {{2, 3}});
    for (int i = 0; i < 9; ++i) {
        raw.points.source.push_back({1, 10 * i});
        raw.points.weight.push_back(static_cast<float>(i));
    }
    raw.points.triangle = {0, 1, 2, 3, 4, 5};
    raw.points.neighbours[0] = {1, 8, 7, -1, 0, 0, 0, 0, 2};
    raw.points.onBorderEdge.resize(9);
    raw.points.onBorderEdge.set(8);
    raw.originalToOutput = {8, 4, 2};
    raw.curveEdges.resize(4, true);
    extract::FinalCleanup::run(raw);
    EXPECT(faceStartOf(raw.mesh), (std::vector<int>{0, 4, 7}));
    EXPECT(raw.mesh.topology().cornerVertex(), (std::vector<int>{0, 1, 2, 3, 3, 2, 4}));
    EXPECT(raw.mesh.pointCount(), 5);
    EXPECT(raw.mesh.points()[4] == (Vec3f{1, 2, 0}), true);
    const auto* kept = raw.mesh.maps().findAs<TypedMap<Color4ub>>(map_names::colors, ElementKind::Point);
    EXPECT(kept != nullptr && kept->size() == 5, true);
    if (kept != nullptr && kept->size() == 5) {
        EXPECT(kept->values()[3].r, std::uint8_t{30});
        EXPECT(kept->values()[4].r, std::uint8_t{80});
    }
    EXPECT(isSharp(raw.mesh, 2, 3), true);
    EXPECT(isSharp(raw.mesh, 0, 1), false);
    EXPECT(raw.points.source.size(), std::size_t{5});
    EXPECT(raw.points.source[4] == (extract::PointSource{1, 80}), true);
    EXPECT(raw.points.weight, (std::vector<float>{0, 1, 2, 3, 8}));
    EXPECT(raw.points.triangle, (std::vector<int>{0, 1, 2, 3}));
    EXPECT(raw.points.neighbours[0], (std::vector<int>{1, 4, -1, -1, 2}));
    EXPECT(raw.points.onBorderEdge.size(), std::size_t{5});
    EXPECT(raw.points.onBorderEdge.test(4), true);
    EXPECT(raw.originalToOutput, (std::vector<int>{4, -1, 2}));
    EXPECT(raw.curveEdges.empty(), true);
}

// Nothing to clean: faces, points and maps stay as they are.
void finalCleanupWithoutChange()
{
    extract::RawMesh raw;
    raw.mesh = polygonMesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, {0, 4}, {0, 1, 2, 3}, FaceMode::Quads);
    raw.points.source.assign(4, extract::PointSource{4, 0});
    extract::FinalCleanup::run(raw);
    EXPECT(raw.mesh.topology().mode(), FaceMode::Quads);
    EXPECT(raw.mesh.pointCount(), 4);
    EXPECT(raw.points.source.size(), std::size_t{4});
}

// The kept border points with an output point get its place back bit for
// bit; a point without mark 1 and a point without output point do not.
void frozenBorderReset()
{
    Mesh triangles = polygonMesh({{0.1f, 0.2f, 0.3f}, {1.7f, 0, 0}, {1.3f, 1.1f, 0}, {0.3f, 1.9f, 0.7f}}, {0, 3, 6},
                                 {0, 1, 2, 0, 2, 3}, FaceMode::Triangles);
    Mesh output = polygonMesh({{9, 9, 9}, {8, 8, 8}, {7, 7, 7}}, {0, 3}, {0, 1, 2}, FaceMode::Polygons);
    const std::vector<FreezeMark> marks{{1}, {1}, {2}, {1}};
    extract::FrozenBorderReset::run(output, triangles, marks, {2, -1, 0, 1});
    EXPECT(output.points()[2] == triangles.points()[0], true);
    EXPECT(output.points()[1] == triangles.points()[3], true);
    EXPECT(output.points()[0] == (Vec3f{9, 9, 9}), true);
}

// Colours by source: a point, the middle of an edge, a point inside the
// triangle with weights 1/2, 1/4, 1/4, a point outside it (one channel below
// 0, one above 255) and a point without source.
void vertexColorTransfer()
{
    Mesh triangles = polygonMesh({{0, 0, 0}, {2, 0, 0}, {0, 2, 0}}, {0, 3}, {0, 1, 2}, FaceMode::Triangles);
    triangles.maps().set(std::make_unique<TypedMap<Color4ub>>(
        std::string(map_names::colors), ElementKind::Point,
        std::vector<Color4ub>{{200, 0, 0, 255}, {0, 200, 0, 255}, {0, 0, 250, 255}}));
    const int edge = triangles.adjacency().edgeBetween(0, 1);
    Mesh output = polygonMesh({{0, 0, 0}, {1, 0, 0}, {0.5f, 0.5f, 0}, {4, 0, 0}, {1, 1, 1}}, {0, 5}, {0, 1, 2, 3, 4},
                              FaceMode::Polygons);
    const std::vector<extract::PointSource> sources{{1, 1}, {2, edge}, {4, 0}, {4, 0}, {0, 0}};
    extract::VertexColorTransfer::run(output, sources, triangles);
    const auto* c = output.maps().findAs<TypedMap<Color4ub>>(map_names::colors, ElementKind::Point);
    EXPECT(c != nullptr && c->size() == 5, true);
    if (c != nullptr && c->size() == 5) {
        EXPECT(c->values()[0] == (Color4ub{0, 200, 0, 255}), true);
        EXPECT(c->values()[1] == (Color4ub{100, 100, 0, 255}), true);
        EXPECT(c->values()[2] == (Color4ub{100, 50, 62, 255}), true);
        EXPECT(c->values()[3] == (Color4ub{0, 255, 0, 255}), true);
        EXPECT(c->values()[4] == (Color4ub{0, 0, 0, 0}), true);
    }
    const extract::SourceWeights w = extract::SourceWeights::of(triangles, output.points(), sources);
    EXPECT(w.rows.starts, (std::vector<int>{0, 1, 3, 6, 9, 9}));
    EXPECT(w.weights[3] == 0.5f && w.weights[4] == 0.25f && w.weights[5] == 0.25f, true);
    // An RGB map stays RGB.
    Mesh rgb = polygonMesh({{0, 0, 0}, {2, 0, 0}, {0, 2, 0}}, {0, 3}, {0, 1, 2}, FaceMode::Triangles);
    rgb.maps().set(std::make_unique<TypedMap<Color3ub>>(std::string(map_names::colors), ElementKind::Point,
                                                        std::vector<Color3ub>{{200, 0, 0}, {0, 200, 0}, {0, 0, 250}}));
    extract::VertexColorTransfer::run(output, sources, rgb);
    const auto* c3 = output.maps().findAs<TypedMap<Color3ub>>(map_names::colors, ElementKind::Point);
    EXPECT(c3 != nullptr && c3->size() == 5 && c3->values()[2] == (Color3ub{100, 50, 62}), true);
}

bool near(Vec3f a, Vec3f b, float tolerance)
{
    return std::fabs(a.x - b.x) <= tolerance && std::fabs(a.y - b.y) <= tolerance && std::fabs(a.z - b.z) <= tolerance;
}

// Input of one relaxation step on `m`: wide length limits, no directions, no
// anisotropy, only `movable` moves.
extract::FinalRelaxation::StepInput stepInput(const Mesh& m, std::initializer_list<int> movable)
{
    extract::FinalRelaxation::StepInput in;
    const auto edges = static_cast<std::size_t>(m.adjacency().edgeCount());
    in.edges.direction.assign(edges, Vec3f{});
    in.edges.lengthMin.assign(edges, 0.0f);
    in.edges.lengthMax.assign(edges, 10.0f);
    in.anisotropy.assign(static_cast<std::size_t>(m.pointCount()), 0.0f);
    in.mobility.assign(static_cast<std::size_t>(m.pointCount()), 0.0f);
    for (const int p : movable) {
        in.mobility[static_cast<std::size_t>(p)] = 1.0f;
    }
    in.lines.resize(edges);
    return in;
}

// One step: an inner point with three edges goes 0.4 of the way to the mean
// of its neighbours; a point with four edges goes 0.33 of the way to the mean
// of the two points on its opposite edge pairs, each weighted by the lengths
// ((0, 1)-(2, 1) at 1.197939, (1, 0)-(1, 2) at 1.096118); a border point
// between two border lines goes along the nearer line.
void relaxationStep()
{
    Mesh fan = polygonMesh({{0.3f, 0.3f, 0}, {1, 0, 0}, {0, 1, 0}, {-1, -1, 0}}, {0, 3, 6, 9}, {0, 1, 2, 0, 2, 3, 0, 3, 1},
                           FaceMode::Polygons);
    std::vector<Vec3f> out = extract::FinalRelaxation::step(fan, stepInput(fan, {0}));
    EXPECT(near(out[0], {0.18f, 0.18f, 0}, 1e-6f), true);
    EXPECT(out[1] == fan.points()[1], true);

    std::vector<Vec3f> p;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            p.push_back({static_cast<float>(i), static_cast<float>(j), 0});
        }
    }
    p[4] = {1.2f, 1.1f, 0};
    Mesh grid = polygonMesh(p, {0, 4, 8, 12, 16}, {0, 1, 4, 3, 1, 2, 5, 4, 3, 4, 7, 6, 4, 5, 8, 7}, FaceMode::Quads);
    out = extract::FinalRelaxation::step(grid, stepInput(grid, {4}));
    const float isoX = (1.197939f + 1.0f) * 0.5f;
    const float isoY = (1.0f + 1.096118f) * 0.5f;
    EXPECT(near(out[4], {isoX * 0.33f + 1.2f * 0.67f, isoY * 0.33f + 1.1f * 0.67f, 0}, 1e-5f), true);

    // Bottom border point (1.3, 0): mean (1, 1/3) laid onto the border line.
    p = {{0, 0, 0}, {1.3f, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}};
    Mesh strip = polygonMesh(p, {0, 4, 8}, {0, 1, 4, 3, 1, 2, 5, 4}, FaceMode::Quads);
    extract::FinalRelaxation::StepInput in = stepInput(strip, {1});
    for (const int e : strip.adjacency().boundaryEdges()) {
        in.lines.set(static_cast<std::size_t>(e));
    }
    out = extract::FinalRelaxation::step(strip, in);
    EXPECT(near(out[1], {1.18f, 0, 0}, 1e-6f), true);
    EXPECT(out[1].y == 0.0f, true);
    // With the inner edge on a line as well the point stays (three lines meet).
    in.lines.set(static_cast<std::size_t>(strip.adjacency().edgeBetween(1, 4)));
    out = extract::FinalRelaxation::step(strip, in);
    EXPECT(out[1] == strip.points()[1], true);
}

// Two quads folded by a right angle along their shared edge.
void relaxationDihedral()
{
    Mesh m = polygonMesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 0, 1}, {1, 1, 1}}, {0, 4, 8},
                         {0, 1, 2, 3, 1, 4, 5, 2}, FaceMode::Quads);
    const std::vector<float> angles = extract::FinalRelaxation::dihedralAngles(m);
    const int shared = m.adjacency().edgeBetween(1, 2);
    EXPECT(std::fabs(angles[static_cast<std::size_t>(shared)] - 1.5707964f) < 1e-6f, true);
    EXPECT(angles[static_cast<std::size_t>(m.adjacency().edgeBetween(0, 1))], 0.0f);
    // A face without a normal (NaN point) gives a NaN angle, not 0.
    Mesh broken = polygonMesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {std::nanf(""), 0, 1}, {1, 1, 1}}, {0, 4, 8},
                              {0, 1, 2, 3, 1, 4, 5, 2}, FaceMode::Quads);
    EXPECT(std::isnan(extract::FinalRelaxation::dihedralAngles(broken)[static_cast<std::size_t>(shared)]), true);
}

// The walk over the square of two triangles (0,0) (2,0) (2,2) (0,2): across
// the diagonal into the other triangle, to a border it cannot cross, along a
// fixed border edge (the shift only counts along it), and from a point.
void surfaceWalk()
{
    Mesh square = polygonMesh({{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}}, {0, 3, 6}, {0, 1, 2, 0, 2, 3},
                              FaceMode::Triangles);
    BitArray none;
    extract::SurfaceWalker free(square, none);
    extract::PointSource s{4, 0};
    Vec3f at{1.5f, 0.5f, 0};
    free.walk(s, at, {-1, 0.5f, 0.3f});
    EXPECT(near(at, {0.5f, 1, 0}, 1e-6f), true);
    EXPECT(s == (extract::PointSource{4, 1}), true);

    s = {4, 0};
    at = {1.5f, 0.5f, 0};
    free.walk(s, at, {2, 0, 0});
    EXPECT(near(at, {2, 0.5f, 0}, 1e-6f), true);
    EXPECT(s == (extract::PointSource{4, 0}), true);

    BitArray border(static_cast<std::size_t>(square.adjacency().edgeCount()));
    for (const int e : square.adjacency().boundaryEdges()) {
        border.set(static_cast<std::size_t>(e));
    }
    extract::SurfaceWalker fixed(square, border);
    const int bottom = square.adjacency().edgeBetween(0, 1);
    s = {2, bottom};
    at = {1, 0, 0};
    fixed.walk(s, at, {0.5f, 0.3f, 0});
    EXPECT(near(at, {1.5f, 0, 0}, 1e-6f), true);
    EXPECT(at.y, 0.0f);
    EXPECT(s == (extract::PointSource{2, bottom}), true);

    // From the corner point (2, 0) into the square: out through the far side of its triangle.
    s = {1, 1};
    at = {2, 0, 0};
    free.walk(s, at, {-0.5f, 0.25f, 0});
    EXPECT(near(at, {1.5f, 0.25f, 0}, 1e-6f), true);
    EXPECT(s.type, 1);
    // A shift without length leaves everything as it is.
    free.walk(s, at, {0, 0, 0});
    EXPECT(near(at, {1.5f, 0.25f, 0}, 1e-6f), true);

    // Floor and wall at a right angle: the walk meets the wall almost head on,
    // the direction laid into the wall nearly vanishes, so it is turned about
    // the shared edge and goes up the wall by the rest of the shift.
    Mesh fold = polygonMesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 0, 1}}, {0, 3, 6}, {0, 1, 2, 2, 1, 3},
                            FaceMode::Triangles);
    extract::SurfaceWalker overFold(fold, none);
    s = {4, 0};
    at = {0.5f, 0.2f, 0};
    overFold.walk(s, at, {1, 0.0002f, 0});
    EXPECT(near(at, {1, 0.2002f, 0.5f}, 1e-5f), true);
    EXPECT(s == (extract::PointSource{4, 1}), true);
}

// A traced square with uniform fields (size 1, combed directions along x and
// y): every edge gets the direction from its first to its second end point
// (or none) and length limits between 0.5 .. 0.7 and 1.5 .. 2.5. The quads
// already have these lengths, so the finished mesh keeps every point (within
// 1e-5, in the plane; the four corners bit for bit, fixed by the bend of the
// border), and nothing is cleaned.
void relaxationOnSquare()
{
    parts::Part part = tracedPart(squareLattice(), {}, [](const Adjacency&, int) { return 0; }, {}, {0});
    part.quadBudget = 9;
    part.edgeLength = 1.0f;
    rawMeshOf(part);
    EXPECT(part.raw.has_value(), true);
    if (!part.raw) {
        return;
    }
    extract::Extraction::quads(part, {});
    const auto triangles = static_cast<std::size_t>(part.mesh.faceCount());
    field::SizeField sizes;
    sizes.sizeU.assign(triangles, 1.0f);
    sizes.sizeV.assign(triangles, 1.0f);
    part.sizeField = std::move(sizes);
    field::CrossField cross;
    cross.frameE1.assign(triangles, Vec3f{1, 0, 0});
    cross.frameE2.assign(triangles, Vec3f{0, 1, 0});
    part.crossField = std::move(cross);
    part.cutGraph->combedU.assign(triangles, {1.0f, 0.0f});
    part.cutGraph->combedV.assign(triangles, {0.0f, 1.0f});

    const extract::FinalRelaxation::EdgeTargets targets = extract::FinalRelaxation::edgeTargets(part);
    const Mesh& out = part.raw->mesh;
    int directed = 0;
    int wrong = 0;
    for (int e = 0; e < out.adjacency().edgeCount(); ++e) {
        const auto i = static_cast<std::size_t>(e);
        const auto ends = out.adjacency().edgeVertices(e);
        Vec3f along = out.points()[static_cast<std::size_t>(ends[1])] - out.points()[static_cast<std::size_t>(ends[0])];
        normalizeInPlace(along);
        const Vec3f r = targets.direction[i];
        if (!(r == Vec3f{})) {
            ++directed;
            wrong += near(r, along, 1e-6f) ? 0 : 1;
        }
        wrong += targets.lengthMin[i] >= 0.5f && targets.lengthMin[i] <= 0.7f ? 0 : 1;
        wrong += targets.lengthMax[i] >= 1.5f && targets.lengthMax[i] <= 2.5f ? 0 : 1;
    }
    EXPECT(wrong, 0);
    EXPECT(directed > 0, true);

    const std::vector<Vec3f> before = out.points();
    std::atomic<bool> cancel{false};
    run::Progress progress(nullptr, nullptr, cancel);
    extract::Extraction::finish(part, {}, progress);
    const std::vector<Vec3f>& after = part.raw->mesh.points();
    EXPECT(after.size(), before.size());
    EXPECT(part.raw->mesh.faceCount(), 9);
    if (after.size() != before.size()) {
        return;
    }
    int moved = 0;
    int offPlane = 0;
    int cornersMoved = 0;
    for (std::size_t i = 0; i < after.size(); ++i) {
        moved += near(after[i], before[i], 1e-5f) ? 0 : 1;
        offPlane += after[i].z == 0.0f ? 0 : 1;
        const bool corner = (before[i].x < 0.01f || before[i].x > 2.99f) && (before[i].y < 0.01f || before[i].y > 2.99f);
        cornersMoved += corner && !(after[i] == before[i]) ? 1 : 0;
    }
    EXPECT(moved, 0);
    EXPECT(offPlane, 0);
    EXPECT(cornersMoved, 0);
}

} // namespace

// Features of the triangle mesh on the quad mesh: face groups without crease
// curves, a crease along the diagonal of a quad (cut, sharp edge, groups per
// face) and a crease across two quads (cut with a new point on the crease
// chain and its source); a closed crease loop with its closing edge, the
// sharp edge decided by the final best candidate (and a sharp edge within
// 6 degrees winning at once), open faces taking the group of their
// neighbours (first one on ties, not across curve edges).
namespace {

// A part with the triangle mesh (points, triangles, face groups), one sharp
// crease along `chain`, and as quad mesh the polygons on the first
// `quadPoints` points, each output point coming from the same point.
parts::Part featurePart(const std::vector<Vec3f>& points, std::vector<int> triangles, std::vector<int> groups,
                        std::vector<int> chain, int quadPoints, std::vector<int> faceStart, std::vector<int> corners)
{
    parts::Part part;
    std::vector<Vec3f> own = points;
    part.mesh.swapPoints(own);
    part.mesh.editTopology().swapCorners(triangles, FaceMode::Triangles);
    part.mesh.topologyChanged();
    if (!groups.empty()) {
        part.mesh.maps().set(
            std::make_unique<TypedMap<int>>(std::string(map_names::faceGroups), ElementKind::Face, std::move(groups)));
    }
    if (!chain.empty()) {
        curves::GuideCurve curve;
        curve.kind = ARM_CURVE_CREASE;
        curve.lineFlags = 1;
        curve.cutVertices = std::move(chain);
        part.curves.push_back(std::move(curve));
    }
    part.trace.emplace(part.mesh.faceCount());
    part.trace->curves = extract::CurveMarks::mark(part);
    extract::RawMesh raw;
    raw.mesh = polygonMesh({points.begin(), points.begin() + quadPoints}, std::move(faceStart), std::move(corners),
                           FaceMode::Polygons);
    const auto count = static_cast<std::size_t>(quadPoints);
    for (std::size_t i = 0; i < count; ++i) {
        raw.points.source.push_back({1, static_cast<int>(i)});
    }
    raw.points.weight.assign(count, 0.0f);
    raw.points.triangle.assign(count, -1);
    for (std::vector<int>& list : raw.points.neighbours) {
        list.assign(count, -1);
    }
    part.raw = std::move(raw);
    return part;
}

const std::vector<Vec3f> featureStrip{{0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}, {2.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f},
                                      {1.0f, 1.0f, 0.0f}, {2.0f, 1.0f, 0.0f}, {1.0f, 0.5f, 0.0f}};

// Without crease curves every quad takes the group of the first triangle.
void featuresWithoutCurves()
{
    parts::Part part = featurePart(featureStrip, {0, 1, 4, 0, 4, 3}, {5, 6}, {}, 6, {0, 4, 8}, {0, 1, 4, 3, 1, 2, 5, 4});
    extract::Extraction::features(part, {});
    const extract::RawMesh& raw = *part.raw;
    EXPECT(faceGroupValues(raw.mesh), (std::vector<int>{5, 5}));
    EXPECT(SharpEdges::edgeBits(raw.mesh).has_value(), false);
    EXPECT(raw.points.onCurve.size(), std::size_t{6});
    EXPECT(raw.points.onCurve.any(), false);
}

// A sharp crease along the diagonal 0-4 of the left quad: the quad is split
// into two triangles along it, the diagonal becomes the one curve edge and
// sharp, the right quad keeps the group of its triangles and comes first.
void featureCreaseAcrossQuad()
{
    parts::Part part = featurePart(featureStrip, {0, 1, 4, 0, 4, 3, 1, 2, 5, 1, 5, 4}, {7, 7, 9, 9}, {0, 4}, 6,
                                   {0, 4, 8}, {0, 1, 4, 3, 1, 2, 5, 4});
    extract::Extraction::features(part, {});
    const extract::RawMesh& raw = *part.raw;
    EXPECT(raw.mesh.pointCount(), 6);
    EXPECT(cornerCounts(raw.mesh), (std::vector<int>{0, 0, 0, 2, 1, 0, 0, 0}));
    EXPECT(raw.points.onCurve.count(), std::size_t{2});
    const Adjacency& adjacency = raw.mesh.adjacency();
    const int diagonal = adjacency.edgeBetween(0, 4);
    EXPECT(diagonal >= 0, true);
    EXPECT(raw.curveEdges.count(), std::size_t{1});
    EXPECT(raw.curveEdges.testBounded(static_cast<std::size_t>(diagonal)), true);
    const std::optional<BitArray> sharp = SharpEdges::edgeBits(raw.mesh);
    EXPECT(sharp.has_value(), true);
    if (sharp) {
        EXPECT(sharp->count(), std::size_t{1});
        EXPECT(sharp->testBounded(static_cast<std::size_t>(diagonal)), true);
    }
    EXPECT(faceGroupValues(raw.mesh), (std::vector<int>{9, 7, 7}));
}

// A sharp crease 0 - 6 - 5 through point 6 = (1, 0.5) of the triangle mesh,
// which the quad mesh does not have: the cut crosses the shared side 1-4 at
// (1, 0.5), takes the new point onto the crease chain there (source: a chain
// edge) and splits both quads into a quad and a triangle.
void featureCreaseThroughQuads()
{
    parts::Part part = featurePart(featureStrip, {0, 1, 6, 0, 6, 3, 3, 6, 4, 1, 2, 6, 6, 2, 5, 6, 5, 4},
                                   {7, 7, 7, 9, 9, 9}, {0, 6, 5}, 6, {0, 4, 8}, {0, 1, 4, 3, 1, 2, 5, 4});
    extract::Extraction::features(part, {});
    const extract::RawMesh& raw = *part.raw;
    EXPECT(raw.mesh.pointCount(), 7);
    EXPECT(cornerCounts(raw.mesh), (std::vector<int>{0, 0, 0, 2, 2, 0, 0, 0}));
    if (raw.mesh.pointCount() != 7) {
        return;
    }
    EXPECT(distance(raw.mesh.points()[6], Vec3f{1.0f, 0.5f, 0.0f}) < 1e-5f, true);
    const extract::PointSource source = raw.points.source[6];
    const Adjacency& triangles = part.mesh.adjacency();
    EXPECT(source.type, 2);
    EXPECT(source.id == triangles.edgeBetween(0, 6) || source.id == triangles.edgeBetween(6, 5), true);
    EXPECT(raw.points.onCurve.testBounded(6), true);
    const Adjacency& adjacency = raw.mesh.adjacency();
    EXPECT(raw.curveEdges.count(), std::size_t{2});
    EXPECT(raw.curveEdges.testBounded(static_cast<std::size_t>(adjacency.edgeBetween(0, 6))), true);
    EXPECT(raw.curveEdges.testBounded(static_cast<std::size_t>(adjacency.edgeBetween(6, 5))), true);
    const std::optional<BitArray> sharp = SharpEdges::edgeBits(raw.mesh);
    EXPECT(sharp.has_value() && sharp->count() == 2, true);
}

// A sharp crease around the middle square 5-6-10-9 of a 3 x 3 quad grid:
// every edge of the closed loop, also the closing one 9-5, becomes a curve
// edge and sharp; no point is added.
void featureCreaseClosedLoop()
{
    std::vector<Vec3f> points;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            points.push_back({static_cast<float>(i), static_cast<float>(j), 0.0f});
        }
    }
    std::vector<int> triangles;
    std::vector<int> faceStart{0};
    std::vector<int> corners;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            const int a = j * 4 + i;
            triangles.insert(triangles.end(), {a, a + 1, a + 5, a, a + 5, a + 4});
            corners.insert(corners.end(), {a, a + 1, a + 5, a + 4});
            faceStart.push_back(static_cast<int>(corners.size()));
        }
    }
    parts::Part part = featurePart(points, triangles, {}, {5, 6, 10, 9, 5}, 16, faceStart, corners);
    extract::Extraction::features(part, {});
    const extract::RawMesh& raw = *part.raw;
    EXPECT(raw.mesh.pointCount(), 16);
    EXPECT(raw.points.onCurve.count(), std::size_t{4});
    const Adjacency& adjacency = raw.mesh.adjacency();
    EXPECT(raw.curveEdges.count(), std::size_t{4});
    const std::optional<BitArray> sharp = SharpEdges::edgeBits(raw.mesh);
    EXPECT(sharp.has_value() && sharp->count() == 4, true);
    for (const auto& [a, b] : {std::pair{5, 6}, std::pair{6, 10}, std::pair{10, 9}, std::pair{9, 5}}) {
        const int e = adjacency.edgeBetween(a, b);
        EXPECT(e >= 0 && raw.curveEdges.testBounded(static_cast<std::size_t>(e)), true);
        EXPECT(e >= 0 && sharp && sharp->testBounded(static_cast<std::size_t>(e)), true);
    }
}

// Whether the quad edge over the triangle points a = 0 and b = 1 = (10, 0, 0)
// becomes sharp. The edges at a are, in point -> edge order, 0-2 at
// `firstDegrees` above and 0-3 at `secondDegrees` below the line a-b; the
// sharp line runs along 0-`sharpEnd`.
bool sharpBetween(float firstDegrees, float secondDegrees, int sharpEnd)
{
    auto towards = [](float degrees, float side) {
        const float r = degrees * 0.017453292f;
        return Vec3f{std::cos(r), side * std::sin(r), 0.0f};
    };
    const Mesh triangles = polygonMesh({{0.0f, 0.0f, 0.0f}, {10.0f, 0.0f, 0.0f}, towards(firstDegrees, 1.0f),
                                        towards(secondDegrees, -1.0f)},
                                       {0, 3, 6}, {0, 2, 3, 3, 2, 1}, FaceMode::Triangles);
    curves::GuideCurve line;
    line.kind = ARM_CURVE_CREASE;
    line.lineFlags = 1;
    line.cutVertices = {0, sharpEnd};
    const std::vector<curves::GuideCurve> lines{line};
    extract::RawMesh raw;
    raw.mesh = polygonMesh({{0.0f, 0.0f, 0.0f}, {10.0f, 0.0f, 0.0f}, {5.0f, 5.0f, 0.0f}}, {0, 3}, {0, 1, 2},
                           FaceMode::Polygons);
    raw.points.source = {{1, 0}, {1, 1}, {0, -1}};
    const int q = raw.mesh.adjacency().edgeBetween(0, 1);
    raw.curveEdges = BitArray(static_cast<std::size_t>(raw.mesh.adjacency().edgeCount()));
    raw.curveEdges.set(static_cast<std::size_t>(q));
    extract::HardEdgeTransfer::run(triangles, lines, raw);
    const std::optional<BitArray> sharp = SharpEdges::edgeBits(raw.mesh);
    return sharp.has_value() && sharp->testBounded(static_cast<std::size_t>(q));
}

// The sharp edge at 8 degrees (counts 4.8) comes first, a better plain edge
// at 4 degrees follows: the final best one decides, not sharp. Against a
// plain edge at 6 degrees the sharp one stays best. A sharp edge within
// 6 degrees wins at once, also after a plain edge at 2 degrees.
void sharpEdgeOfBestCandidate()
{
    EXPECT(sharpBetween(8.0f, 4.0f, 2), false);
    EXPECT(sharpBetween(8.0f, 6.0f, 2), true);
    EXPECT(sharpBetween(2.0f, 5.0f, 3), true);
}

// Three quads in a row: the outer ones take the groups 3 and 8 of the
// triangles their corners lie in, the middle one has corners in both and
// stays open. Its neighbours tie (8 across 2-6 comes first, then 3 across
// 5-1): 8. With 2-6 a curve edge only 3 counts.
void openFacesTakeNeighbourGroup()
{
    Mesh triangles = polygonMesh({{0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {1.0f, 1.0f, 0.0f}},
                                 {0, 3, 6}, {0, 1, 2, 2, 1, 3}, FaceMode::Triangles);
    triangles.maps().set(std::make_unique<TypedMap<int>>(std::string(map_names::faceGroups), ElementKind::Face,
                                                         std::vector<int>{3, 8}));
    auto run = [&](bool curveAcross) {
        extract::RawMesh raw;
        std::vector<Vec3f> points;
        for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 4; ++i) {
                points.push_back({static_cast<float>(i), static_cast<float>(j), 0.0f});
            }
        }
        raw.mesh = polygonMesh(points, {0, 4, 8, 12}, {0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6}, FaceMode::Polygons);
        for (const int v : {0, 1, 2, 3, 4, 5, 6, 7}) {
            raw.points.source.push_back({4, (v % 4) < 2 ? 0 : 1});
        }
        const Adjacency& adjacency = raw.mesh.adjacency();
        raw.curveEdges = BitArray(static_cast<std::size_t>(adjacency.edgeCount()));
        if (curveAcross) {
            raw.curveEdges.set(static_cast<std::size_t>(adjacency.edgeBetween(2, 6)));
        }
        extract::GroupTransfer::run(triangles, raw);
        return faceGroupValues(raw.mesh);
    };
    EXPECT(run(false), (std::vector<int>{3, 8, 8}));
    EXPECT(run(true), (std::vector<int>{3, 3, 8}));
}

} // namespace

// The result: host groups back, points split along sharp edges, parts put
// together with their maps and new face groups, the mirror at the end.
namespace {

void setFaceInts(Mesh& m, std::string_view name, std::vector<int> values)
{
    m.maps().set(std::make_unique<TypedMap<int>>(std::string(name), ElementKind::Face, std::move(values)));
}

const TypedMap<int>* faceIntsOf(const Mesh& m, std::string_view name)
{
    return m.maps().findAs<TypedMap<int>>(name, ElementKind::Face);
}

// Three quads in a row with the regions 0, 2 and 5 (5 lies beyond the table).
Mesh threeRegions(std::vector<int> regions)
{
    Mesh m = polygonMesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}, {3, 1, 0}},
                         {0, 4, 8, 12}, {0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6}, FaceMode::Polygons);
    setFaceInts(m, map_names::faceGroups, std::move(regions));
    return m;
}

// Material, face group and smoothing groups by the region a face had; the
// smoothing groups do not take the host face group as region; a region
// beyond the table gives 0; materials that were the groups come back as the
// regions; a map shorter than the faces changes nothing.
void groupRestoreRules()
{
    prep::GroupRestoreTable table;
    table.groupMaterial = {10, 11, 12};
    table.groupFaceGroup = {1, 2, 0};
    table.groupSmoothing = {1u, 2u, 4u};
    table.restore = true;
    Mesh m = threeRegions({0, 2, 5});
    merge::GroupRestore::run(m, table);
    EXPECT(faceIntsOf(m, map_names::materials) != nullptr, true);
    if (faceIntsOf(m, map_names::materials) != nullptr) {
        EXPECT(faceIntsOf(m, map_names::materials)->values(), (std::vector<int>{10, 12, 0}));
    }
    EXPECT(faceGroupValues(m), (std::vector<int>{1, 0, 0}));
    const auto* smoothing = m.maps().findAs<TypedMap<std::uint32_t>>(map_names::smoothingGroups, ElementKind::Face);
    EXPECT((smoothing != nullptr && smoothing->values() == std::vector<std::uint32_t>{1u, 4u, 0u}), true);

    prep::GroupRestoreTable copies;
    copies.materialsAreGroups = true;
    copies.restore = true;
    Mesh n = threeRegions({0, 2, 5});
    merge::GroupRestore::run(n, copies);
    EXPECT(faceIntsOf(n, map_names::materials) != nullptr, true);
    if (faceIntsOf(n, map_names::materials) != nullptr) {
        EXPECT(faceIntsOf(n, map_names::materials)->values(), (std::vector<int>{0, 2, 5}));
    }
    EXPECT(faceGroupValues(n), (std::vector<int>{0, 2, 5}));

    Mesh shortMap = threeRegions({0, 2});
    merge::GroupRestore::run(shortMap, table);
    EXPECT(faceIntsOf(shortMap, map_names::materials) == nullptr, true);
    EXPECT(faceGroupValues(shortMap), (std::vector<int>{0, 2}));
}

// Six quads on 4 x 3 points (row y = 0: points 0..3, y = 1: 4..7, y = 2:
// 8..11) with a point map holding the point number.
Mesh sixQuads()
{
    std::vector<Vec3f> points;
    for (int y = 0; y < 3; ++y) {
        for (int x = 0; x < 4; ++x) {
            points.push_back({static_cast<float>(x), static_cast<float>(y), 0.0f});
        }
    }
    Mesh m = polygonMesh(points, {0, 4, 8, 12, 16, 20, 24},
                         {0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 4, 5, 9, 8, 5, 6, 10, 9, 6, 7, 11, 10}, FaceMode::Polygons);
    std::vector<float> numbers(12);
    std::iota(numbers.begin(), numbers.end(), 0.0f);
    m.maps().set(std::make_unique<TypedMap<float>>(std::string(map_names::density), ElementKind::Point, numbers));
    return m;
}

// A sharp line 1 - 5 - 9 across the grid: the candidates are 5 (two sharp
// edges), then 1 and 9 (a sharp edge and the border); each has two fans,
// the second one (with quads 1 and 4) gets the new points 12 (for 5), 13 (for
// 1) and 14 (for 9). Both sides of the line stay sharp. A line 1 - 5 that
// ends inside splits only its border point.
void creaseSplitFans()
{
    Mesh m = sixQuads();
    setSharp(m, {{1, 5}, {5, 9}});
    EXPECT(merge::CreaseSplit::run(m), 3);
    EXPECT(m.pointCount(), 15);
    EXPECT(m.topology().cornerVertex(),
           (std::vector<int>{0, 1, 5, 4, 13, 2, 6, 12, 2, 3, 7, 6, 4, 5, 9, 8, 12, 6, 10, 14, 6, 7, 11, 10}));
    EXPECT(m.points()[12] == m.points()[5] && m.points()[13] == m.points()[1] && m.points()[14] == m.points()[9], true);
    const auto* numbers = m.maps().findAs<TypedMap<float>>(map_names::density, ElementKind::Point);
    EXPECT(numbers != nullptr && numbers->values().size() == 15, true);
    if (numbers != nullptr && numbers->values().size() == 15) {
        EXPECT(numbers->values()[12] == 5.0f && numbers->values()[13] == 1.0f && numbers->values()[14] == 9.0f, true);
    }
    EXPECT(isSharp(m, 1, 5) && isSharp(m, 5, 9) && isSharp(m, 13, 12) && isSharp(m, 12, 14), true);
    EXPECT(isSharp(m, 12, 6), false);

    Mesh inside = sixQuads();
    setSharp(inside, {{1, 5}});
    EXPECT(merge::CreaseSplit::run(inside), 1);
    EXPECT(v(inside.topology().cornerVertex()).at(4), 12);
    EXPECT(v(inside.topology().cornerVertex()).at(7), 5);

    Mesh smooth = sixQuads();
    EXPECT(merge::CreaseSplit::run(smooth), 0);
}

parts::Part finishedPart(int index, Mesh mesh)
{
    parts::Part part;
    part.index = index;
    extract::RawMesh raw;
    raw.mesh = std::move(mesh);
    part.raw = std::move(raw);
    return part;
}

// Three parts: a quad of group 0, a skipped part, and two quads whose face
// group map is empty and which have the only sharp edge. The empty map gets
// group 1 (0 is taken), the sharp edge comes over with the right number; a
// half edge map of the first part only grows by the corners of the third.
// A first part without points and faces hands over to the next one.
void resultMergeMaps()
{
    std::vector<parts::Part> list;
    Mesh first = polygonMesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, {0, 4}, {0, 1, 2, 3}, FaceMode::Polygons);
    setFaceInts(first, map_names::faceGroups, {0});
    first.maps().set(std::make_unique<TypedMap<int>>("marks", ElementKind::HalfEdge, std::vector<int>{1, 2, 3, 4}));
    list.push_back(finishedPart(0, std::move(first)));
    parts::Part skipped;
    skipped.index = 1;
    skipped.skipped = true;
    list.push_back(std::move(skipped));
    Mesh third = strip();
    setFaceInts(third, map_names::faceGroups, {});
    setSharp(third, {{1, 4}});
    list.push_back(finishedPart(2, std::move(third)));
    prep::Prepared prepared;
    Mesh result = merge::ResultAssembly::run(list, prepared, {});
    EXPECT(result.pointCount(), 10);
    EXPECT(result.faceCount(), 3);
    EXPECT(result.topology().mode(), FaceMode::Polygons);
    EXPECT(faceGroupValues(result), (std::vector<int>{0, 1, 1}));
    const std::optional<BitArray> sharp = SharpEdges::edgeBits(result);
    EXPECT(sharp.has_value() && sharp->size() == static_cast<std::size_t>(result.adjacency().edgeCount()), true);
    EXPECT(isSharp(result, 5, 8), true);
    EXPECT(sharp.has_value() && sharp->count() == 1, true);
    const auto* marks = result.maps().findAs<TypedMap<int>>("marks", ElementKind::HalfEdge);
    EXPECT(marks != nullptr ? marks->values() : std::vector<int>{},
           (std::vector<int>{1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0}));

    std::vector<parts::Part> handOver;
    Mesh empty;
    setFaceInts(empty, map_names::materials, {});
    handOver.push_back(finishedPart(0, std::move(empty)));
    handOver.push_back(finishedPart(1, strip()));
    Mesh second = merge::ResultAssembly::run(handOver, prepared, {});
    EXPECT(second.faceCount(), 2);
    EXPECT(faceIntsOf(second, map_names::materials) == nullptr, true);

    std::vector<parts::Part> none;
    none.push_back(finishedPart(0, Mesh{}));
    bool failed = false;
    try {
        (void)merge::ResultAssembly::run(none, prepared, {});
    } catch (const api::StatusError& error) {
        failed = error.status() == ARM_ERROR_NO_RESULT;
    }
    EXPECT(failed, true);
}

// Mirror at the end on the plane x = 0 of the preparation (the middle of the
// box would be x = -1): the points on the plane are shared, the copies of the
// others follow in point order, the mirror faces run the other way, face
// groups and the sharp edge follow.
void resultMirror()
{
    Mesh strip = polygonMesh({{-2, 0, 0}, {-1, 0, 0}, {0, 0, 0}, {-2, 1, 0}, {-1, 1, 0}, {0, 1, 0}}, {0, 4, 8},
                             {0, 1, 4, 3, 1, 2, 5, 4}, FaceMode::Polygons);
    setFaceInts(strip, map_names::faceGroups, {2, 3});
    setSharp(strip, {{1, 4}});
    std::vector<parts::Part> list;
    list.push_back(finishedPart(0, std::move(strip)));
    prep::Prepared prepared;
    prepared.state.symmetry.on = true;
    prepared.state.symmetry.mirrorAtEnd = true;
    prepared.state.symmetry.axis[0] = true;
    prepared.state.symmetryPlane = {0.0f, 0.0f, 0.0f};
    Mesh result = merge::ResultAssembly::run(list, prepared, {});
    EXPECT(result.pointCount(), 10);
    EXPECT(result.topology().cornerVertex(), (std::vector<int>{0, 1, 4, 3, 1, 2, 5, 4, 8, 9, 7, 6, 9, 5, 2, 7}));
    EXPECT(result.points().size() == 10 && result.points()[6] == (Vec3f{2, 0, 0}) && result.points()[9] == (Vec3f{1, 1, 0}),
           true);
    EXPECT(faceGroupValues(result), (std::vector<int>{2, 3, 2, 3}));
    EXPECT(isSharp(result, 1, 4) && isSharp(result, 7, 9), true);
    EXPECT(SharpEdges::edgeBits(result).has_value() && SharpEdges::edgeBits(result)->count() == 2, true);
}

// Points and corners of the result of three separate cubes of different
// size with these settings (the second cube a little moved with `moved`), on
// a context using `cache` when given; `reused` receives whether the run took
// its preparation from the cache.
std::pair<std::vector<float>, std::vector<int>> threeCubesWith(const arm_settings& settings, arm_prep_cache* cache,
                                                               int* reused = nullptr, bool moved = false)
{
    std::vector<float> points;
    std::vector<int> sizes;
    std::vector<int> corners;
    const int cube[6][4] = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
    for (int k = 0; k < 3; ++k) {
        const float side = 1.0f + 0.5f * static_cast<float>(k);
        const float x = 4.0f * static_cast<float>(k) + (moved && k == 1 ? 0.25f : 0.0f);
        for (int p = 0; p < 8; ++p) {
            points.insert(points.end(), {x + ((p == 1 || p == 2 || p == 5 || p == 6) ? side : 0.0f),
                                         (p == 2 || p == 3 || p == 6 || p == 7) ? side : 0.0f, p >= 4 ? side : 0.0f});
        }
        for (const auto& face : cube) {
            sizes.push_back(4);
            for (const int c : face) {
                corners.push_back(8 * k + c);
            }
        }
    }
    arm_context* context = nullptr;
    arm_context_create(&context);
    arm_context_use_prep_cache(context, cache);
    arm_settings_set(context, &settings);
    arm_input_mesh(context, static_cast<int32_t>(points.size() / 3), points.data(), static_cast<int32_t>(sizes.size()),
                   sizes.data(), static_cast<int32_t>(corners.size()), corners.data());
    std::pair<std::vector<float>, std::vector<int>> out;
    arm_stage stage = ARM_STAGE_INPUT;
    arm_mesh_view view{};
    if (arm_run(context, ARM_STAGE_RESULT, nullptr, nullptr) == ARM_OK &&
        arm_stage_mesh(context, &stage, &view) == ARM_OK) {
        out.first.assign(view.points_xyz, view.points_xyz + 3 * static_cast<std::ptrdiff_t>(view.point_count));
        out.second.assign(view.corner_points, view.corner_points + view.corner_count);
    }
    if (reused != nullptr) {
        int32_t flag = -1;
        *reused = arm_stage_preparation_reused(context, &flag) == ARM_OK ? flag : -1;
    }
    arm_context_destroy(context);
    return out;
}

arm_settings cubeSettings(int threads)
{
    arm_settings settings;
    arm_settings_default(&settings);
    settings.quad_count = 600;
    settings.thread_count = threads;
    return settings;
}

// Points and corners of the result of the three cubes, computed with `threads` threads.
std::pair<std::vector<float>, std::vector<int>> threeCubes(int threads)
{
    return threeCubesWith(cubeSettings(threads), nullptr);
}

bool sameResult(const std::pair<std::vector<float>, std::vector<int>>& a,
                const std::pair<std::vector<float>, std::vector<int>>& b)
{
    return !a.second.empty() && a.second == b.second && a.first.size() == b.first.size() &&
           std::memcmp(a.first.data(), b.first.data(), a.first.size() * sizeof(float)) == 0;
}

// The parts on several threads give the result of one thread, bit for bit.
void partsOnThreads()
{
    const auto one = threeCubes(1);
    const auto three = threeCubes(3);
    EXPECT(one.second.empty(), false);
    EXPECT(one.second, three.second);
    EXPECT(one.first.size() == three.first.size() &&
               std::memcmp(one.first.data(), three.first.data(), one.first.size() * sizeof(float)) == 0,
           true);
}

// The remembered preparation: a run with a setting the preparation does not
// read takes it from the cache and gives the result computed without the
// cache, bit for bit; another quad count or a moved point is no hit (and gives
// the computed result as well). The cache outlives the host's hold while a
// context uses it.
void rememberedPreparation()
{
    arm_prep_cache* cache = nullptr;
    EXPECT(arm_prep_cache_create(&cache), ARM_OK);
    arm_settings settings = cubeSettings(3);
    int reused = -1;
    const auto first = threeCubesWith(settings, cache, &reused);
    EXPECT(reused, 0);
    EXPECT(sameResult(first, threeCubesWith(settings, nullptr)), true);

    settings.curvature = 80.0f;
    settings.clean_output = 1;
    const auto curvature = threeCubesWith(settings, cache, &reused);
    EXPECT(reused, 1);
    EXPECT(sameResult(curvature, threeCubesWith(settings, nullptr)), true);

    settings.quad_count = 900;
    const auto count = threeCubesWith(settings, cache, &reused);
    EXPECT(reused, 0);
    EXPECT(sameResult(count, threeCubesWith(settings, nullptr)), true);

    const auto moved = threeCubesWith(settings, cache, &reused, true);
    EXPECT(reused, 0);
    EXPECT(sameResult(moved, threeCubesWith(settings, nullptr, nullptr, true)), true);

    arm_context* context = nullptr;
    arm_context_create(&context);
    EXPECT(arm_context_use_prep_cache(context, cache), ARM_OK);
    arm_prep_cache_destroy(cache);
    EXPECT(arm_context_use_prep_cache(context, nullptr), ARM_OK);
    arm_context_destroy(context);
}

} // namespace

int main()
{
    numberingTriangles();
    numberingQuads();
    numberingPolygons();
    conversions();
    mapMixing();
    edgeChains();
    meshChanges();
    replacementNormals();
    inputBuilder();
    introsortTies();
    costQueueBins();
    groupUnionRegions();
    meshCleaning();
    flatTriangleFlip();
    requiredFlatFlip();
    vertexFan();
    creaseLimitInDegrees();
    creaseChainBreak();
    lineRestore();
    finalLines();
    smootherRepair();
    mirrorCut();
    mirrorDetection();
    mirrorLineBend();
    partsAndBudget();
    curveCut();
    borderLockAxes();
    directSolverOrders();
    driverFixing();
    driverFixingLists();
    driverLocalUpdate();
    jacobiEigen();
    curvatureCylinder();
    crossFieldSphere();
    crossFieldTorus();
    eliminationColumnOrder();
    eliminationRules();
    eliminationColumnBlocks();
    constraintTermRules();
    constraintPeriodRules();
    chainStartRules();
    chainsOnOpenMesh();
    curveChainExclusion();
    curveChainClosedRun();
    systemTransitionTerms();
    systemPenaltyRows();
    systemWithoutUnknown();
    systemOnDisc();
    systemTransitions();
    solutionOnDisc();
    solutionRounds();
    solutionTransitions();
    reweightingFactors();
    mirrorCorners();
    solutionSeeds();
    smoothingNeighbourRules();
    reweightingAveraging();
    stopInRound();
    finalRoundingValues();
    valuesAfterFailedElimination();
    chainEndsFailed();
    mirroredCoordinates();
    borderRoundingSelection();
    borderRoundingRounds();
    bordersBehindRegions();
    borderRoundingSharedPoint();
    bordersWithCurves();
    regionsAndBridges();
    softSystemRows();
    bridgeRowClasses();
    bridgeFieldSizes();
    regionSelection();
    regionRounds();
    regionsAfterFixedPositions();
    seamRowRules();
    cleanupFirstPass();
    cleanupWithOrigins();
    cleanupKeepingMaps();
    doubletAtBorder();
    doubletInside();
    dissolveStraightPoint();
    dissolveRing();
    holeSplitFans();
    holesAndCollapse();
    holeNotches();
    holeGapRings();
    traceSquare();
    traceStar();
    traceCylinder();
    traceFanOtherWay();
    traceWalkLimits();
    traceReseed();
    tracePlaceWithoutNumber();
    ngonCornerFlags();
    ngonRules();
    ngonStripAndGrid();
    ngonNewPointParents();
    ngonParity();
    ngonDoubletsAndSmoothing();
    ngonWholeRun();
    ngonTwoCommits();
    faceClosingRules();
    faceClosingBranches();
    splitIntoQuads();
    borderClassification();
    rawSquare();
    rawStar();
    rawCylinder();
    rawFrozenPlate();
    rawFrozenPlateRings();
    finalCleanupMapsAndUnused();
    finalCleanupWithoutChange();
    frozenBorderReset();
    vertexColorTransfer();
    relaxationStep();
    relaxationDihedral();
    surfaceWalk();
    relaxationOnSquare();
    featuresWithoutCurves();
    featureCreaseAcrossQuad();
    featureCreaseThroughQuads();
    groupRestoreRules();
    creaseSplitFans();
    resultMergeMaps();
    resultMirror();
    featureCreaseClosedLoop();
    sharpEdgeOfBestCandidate();
    openFacesTakeNeighbourGroup();
    partsOnThreads();
    rememberedPreparation();
    std::printf("Pruefungen: %d, Fehler: %d\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
