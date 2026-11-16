configure_file(src/api/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/generated/version.h @ONLY)

set(ARM_FREMD ${PROJECT_SOURCE_DIR}/fremd)

# Polygon tessellator (third-party C sources, unchanged): same floating point
# rules as the core, its own code is not checked for warnings.
add_library(arm_tess OBJECT
    ${ARM_FREMD}/glu-9.0.3/src/libtess/dict.c
    ${ARM_FREMD}/glu-9.0.3/src/libtess/geom.c
    ${ARM_FREMD}/glu-9.0.3/src/libtess/memalloc.c
    ${ARM_FREMD}/glu-9.0.3/src/libtess/mesh.c
    ${ARM_FREMD}/glu-9.0.3/src/libtess/normal.c
    ${ARM_FREMD}/glu-9.0.3/src/libtess/priorityq.c
    ${ARM_FREMD}/glu-9.0.3/src/libtess/render.c
    ${ARM_FREMD}/glu-9.0.3/src/libtess/sweep.c
    ${ARM_FREMD}/glu-9.0.3/src/libtess/tess.c
    ${ARM_FREMD}/glu-9.0.3/src/libtess/tessmono.c)
target_include_directories(arm_tess PRIVATE
    src/tess/gl ${ARM_FREMD}/glu-9.0.3/include ${ARM_FREMD}/glu-9.0.3/src/include)
target_link_libraries(arm_tess PRIVATE arm_rechenregeln)
target_compile_options(arm_tess PRIVATE -w)
set_target_properties(arm_tess PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    C_VISIBILITY_PRESET hidden)

# Sources are compiled once and packed into both library kinds.
add_library(arekoremesher_objekte OBJECT
    src/api/c_api.cpp
    src/api/context.cpp
    src/api/stage_view.cpp
    src/api/status.cpp
    src/curves/closest_face.cpp
    src/curves/curve_distribution.cpp
    src/curves/curve_projector.cpp
    src/curves/polyline_locator.cpp
    src/curves/segment_cutter.cpp
    src/curves/segment_preparation.cpp
    src/extract/border_classifier.cpp
    src/extract/curve_cut_transfer.cpp
    src/extract/curve_marks.cpp
    src/extract/doublets.cpp
    src/extract/edge_dissolve.cpp
    src/extract/extraction.cpp
    src/extract/face_closing.cpp
    src/extract/face_walker.cpp
    src/extract/feature_transfer.cpp
    src/extract/group_transfer.cpp
    src/extract/hard_edge_transfer.cpp
    src/extract/hole_closing.cpp
    src/extract/isoline_tracer.cpp
    src/extract/ngon_corner_classifier.cpp
    src/extract/ngon_cut_check.cpp
    src/extract/ngon_finish.cpp
    src/extract/ngon_new_points.cpp
    src/extract/ngon_optimizer.cpp
    src/extract/ngon_parity_balance.cpp
    src/extract/ngon_polygon_list.cpp
    src/extract/ngon_polygon_splitter.cpp
    src/extract/ngon_split_round.cpp
    src/extract/ngon_split_templates.cpp
    src/extract/ngon_strip_grid.cpp
    src/extract/node_search.cpp
    src/extract/output_points.cpp
    src/extract/quad_mesh_builder.cpp
    src/extract/seed_sectors.cpp
    src/extract/segment_tracer.cpp
    src/extract/trace.cpp
    src/extract/trace_frame.cpp
    src/extract/trace_nodes.cpp
    src/extract/final_cleanup.cpp
    src/extract/final_relaxation.cpp
    src/extract/frozen_border_reset.cpp
    src/extract/source_weights.cpp
    src/extract/surface_walker.cpp
    src/extract/vertex_color_transfer.cpp
    src/field/border_lock.cpp
    src/field/border_sizes.cpp
    src/field/cross_field_rounding.cpp
    src/field/cross_field_solver.cpp
    src/field/cross_field_system.cpp
    src/field/cross_frames.cpp
    src/field/curve_alignment.cpp
    src/field/curve_edges.cpp
    src/field/curvature_field.cpp
    src/field/curvature_radii.cpp
    src/field/curve_influence.cpp
    src/field/cut_copy.cpp
    src/field/cut_graph_builder.cpp
    src/field/density_factors.cpp
    src/field/disk_cut.cpp
    src/field/dual_spanning_tree.cpp
    src/field/edge_curvature.cpp
    src/field/field_combing.cpp
    src/field/part_edge_length.cpp
    src/field/singularity_finder.cpp
    src/field/size_field_builder.cpp
    src/field/size_smoothing.cpp
    src/field/steiner_tree.cpp
    src/field/symmetric_eigen3.cpp
    src/input/input_builder.cpp
    src/input/mesh_input.cpp
    src/mesh/adjacency.cpp
    src/mesh/attribute_map.cpp
    src/mesh/attribute_maps.cpp
    src/mesh/bit_array.cpp
    src/mesh/boundary_loops.cpp
    src/mesh/connected_parts.cpp
    src/mesh/csr.cpp
    src/mesh/degenerate_cleanup.cpp
    src/mesh/edge_chains.cpp
    src/mesh/edge_groups.cpp
    src/mesh/face_rebuild.cpp
    src/mesh/link_path.cpp
    src/mesh/measures.cpp
    src/mesh/mesh.cpp
    src/mesh/modification.cpp
    src/mesh/normals.cpp
    src/mesh/sharp_edges.cpp
    src/mesh/topology.cpp
    src/mesh/vertex_fan.cpp
    src/merge/crease_split.cpp
    src/merge/group_restore.cpp
    src/merge/result_assembly.cpp
    src/param/anchor_constraints.cpp
    src/param/border_rounding.cpp
    src/param/bridge_rows.cpp
    src/param/chain_builder.cpp
    src/param/chain_sides.cpp
    src/param/constraint_periods.cpp
    src/param/constraint_set.cpp
    src/param/constraint_terms.cpp
    src/param/curve_chain_builder.cpp
    src/param/elimination.cpp
    src/param/elimination_columns.cpp
    src/param/energy_rows.cpp
    src/param/final_rounding.cpp
    src/param/iso_line_constraints.cpp
    src/param/region_bridges.cpp
    src/param/region_growth.cpp
    src/param/region_rounding.cpp
    src/param/reweighting.cpp
    src/param/seam_rows.cpp
    src/param/soft_system.cpp
    src/param/solution_builder.cpp
    src/param/solution_seeds.cpp
    src/param/solve_strategy.cpp
    src/param/system_builder.cpp
    src/param/transition_constraints.cpp
    src/param/unknown_numbering.cpp
    src/param/uv_output.cpp
    src/param/variable_flags.cpp
    src/parts/part.cpp
    src/parts/part_builder.cpp
    src/parts/part_curve_cut.cpp
    src/prep/cad_group_merge.cpp
    src/prep/cleanup.cpp
    src/prep/close_point_merge.cpp
    src/prep/cost_queue.cpp
    src/prep/decimator.cpp
    src/prep/degenerate_groups.cpp
    src/prep/denoiser.cpp
    src/prep/duplicate_faces.cpp
    src/prep/edge_dissolver.cpp
    src/prep/edge_splitter.cpp
    src/prep/face_borders.cpp
    src/prep/feature_chains.cpp
    src/prep/feature_lines.cpp
    src/prep/flat_triangle_check.cpp
    src/prep/freeze_marks.cpp
    src/prep/group_union.cpp
    src/prep/hole_filler.cpp
    src/prep/import_attributes.cpp
    src/prep/line_refresh.cpp
    src/prep/mesh_cleaner.cpp
    src/prep/mesh_repair.cpp
    src/prep/net_lines.cpp
    src/prep/normal_split_edges.cpp
    src/prep/part_symmetry.cpp
    src/prep/polygon_triangulator.cpp
    src/prep/preparation.cpp
    src/prep/reference_length.cpp
    src/prep/refiner.cpp
    src/prep/seam_welder.cpp
    src/prep/sliver_fixer.cpp
    src/prep/small_region_merge.cpp
    src/prep/smoother.cpp
    src/prep/symmetry.cpp
    src/prep/symmetry_cache.cpp
    src/prep/symmetry_cut.cpp
    src/prep/symmetry_detector.cpp
    src/prep/triangulator.cpp
    src/run/pipeline.cpp
    src/run/prep_cache.cpp
    src/run/progress.cpp
    src/settings/settings.cpp
    src/solver/active_set_cg.cpp
    src/solver/constrained_driver.cpp
    src/solver/direct_solver.cpp
    src/solver/normal_equations.cpp
    src/solver/sparse_matrix.cpp
    src/util/index_queue.cpp)
target_include_directories(arekoremesher_objekte
    PUBLIC include
    PRIVATE src src/api ${CMAKE_CURRENT_BINARY_DIR}/generated)
target_include_directories(arekoremesher_objekte SYSTEM PRIVATE
    src/tess/gl ${ARM_FREMD}/cdt-1.3.6/include ${ARM_FREMD}/glu-9.0.3/include ${ARM_FREMD}/eigen-3.4.0)
target_link_libraries(arekoremesher_objekte PUBLIC arm_rechenregeln)
set_target_properties(arekoremesher_objekte PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON)

# The parts of a mesh are computed on several threads.
find_package(Threads REQUIRED)

# Shared library for Blender (ctypes).
add_library(arekoremesher SHARED $<TARGET_OBJECTS:arekoremesher_objekte> $<TARGET_OBJECTS:arm_tess>)
target_include_directories(arekoremesher PUBLIC include)
target_link_libraries(arekoremesher PRIVATE Threads::Threads)
target_link_options(arekoremesher PRIVATE -Wl,--no-undefined -Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/exports.map)
set_target_properties(arekoremesher PROPERTIES LINK_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/exports.map)

# Static library (PIC) for hosts that link the core into their own module.
add_library(arekoremesher_statisch STATIC $<TARGET_OBJECTS:arekoremesher_objekte> $<TARGET_OBJECTS:arm_tess>)
target_include_directories(arekoremesher_statisch PUBLIC include)
target_compile_definitions(arekoremesher_statisch INTERFACE ARM_STATIC)
target_link_libraries(arekoremesher_statisch INTERFACE Threads::Threads)
set_target_properties(arekoremesher_statisch PROPERTIES OUTPUT_NAME arekoremesher)

# Check program for the core logic (run by werkzeuge/bauen.sh after each build).
add_executable(arekoremesher_pruefung pruefung/pruefung.cpp $<TARGET_OBJECTS:arekoremesher_objekte> $<TARGET_OBJECTS:arm_tess>)
target_include_directories(arekoremesher_pruefung PRIVATE include src src/api ${CMAKE_CURRENT_BINARY_DIR}/generated)
target_link_libraries(arekoremesher_pruefung PRIVATE arm_rechenregeln Threads::Threads)
