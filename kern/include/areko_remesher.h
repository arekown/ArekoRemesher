/*
 * ArekoRemesher core - public C interface.
 *
 * This is the only public header of the core library. It is plain C so that
 * Blender (Python/ctypes) and other hosts with their own C++ standard library
 * can bind the same binary without sharing a C++ ABI.
 *
 * Rules:
 *  - No exception ever crosses this interface; failures are reported as
 *    arm_status codes plus a message stored in the context.
 *  - There is no global state; everything hangs off an arm_context or an
 *    arm_prep_cache the host owns.
 *  - Memory returned by the core is owned by the core. Pointers returned by
 *    the arm_stage_* queries stay valid until the next arm_run, arm_input_mesh
 *    or arm_context_destroy on the same context.
 *  - All calls on one context come from one thread, except arm_cancel. An
 *    arm_prep_cache may be used by several contexts on several threads.
 */
#ifndef AREKO_REMESHER_H
#define AREKO_REMESHER_H

#include <stdint.h>

#if defined(ARM_STATIC)
#define ARM_API
#elif defined(__GNUC__) || defined(__clang__)
#define ARM_API __attribute__((visibility("default")))
#else
#define ARM_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Result code of every fallible call. Values are stable. */
typedef enum arm_status {
    ARM_OK = 0,
    ARM_ERROR_INVALID_ARGUMENT = 1,
    ARM_ERROR_OUT_OF_MEMORY = 2,
    ARM_ERROR_INTERNAL = 3,
    ARM_ERROR_NO_MESH = 4,             /* no input mesh or a mesh without faces */
    ARM_ERROR_NO_RESULT = 5,           /* the result has no faces */
    ARM_CANCELLED = 6,                 /* arm_cancel stopped the run */
    ARM_ERROR_INVALID_SETTINGS = 7,    /* e.g. neither a quad count nor an edge length */
    ARM_ERROR_STAGE_NOT_AVAILABLE = 8  /* stage not built yet or no stage computed */
} arm_status;

/* Stages a run can stop at (insight views). Values are stable. */
typedef enum arm_stage {
    ARM_STAGE_INPUT = 0,
    ARM_STAGE_PREPARED = 1,
    ARM_STAGE_CLEANED = 2,
    ARM_STAGE_MIRROR_HALF = 3,
    ARM_STAGE_FEATURE_LINES = 4,
    ARM_STAGE_PARTS = 5,
    ARM_STAGE_CURVES = 6,
    ARM_STAGE_CURVATURE = 7,
    ARM_STAGE_CROSS_FIELD = 8,
    ARM_STAGE_CUTS = 9,
    ARM_STAGE_SIZE_FIELD = 10,
    ARM_STAGE_PARAMETERIZATION = 11,
    ARM_STAGE_RAW_MESH = 12,           /* stage mesh: polygons of the grid faces, all parts */
    ARM_STAGE_QUADS = 13,              /* stage mesh: the raw mesh with its larger faces cut into quads */
    ARM_STAGE_RESULT = 14              /* finished quad mesh, sharp edges, colors (keep_colors), no lines or parts;
                                        * per face material and smoothing groups where the input had them, face group
                                        * (input with only materials: material; only smoothing groups: region number) */
} arm_stage;

/* Integer values per face. */
typedef enum arm_face_values {
    ARM_FACE_MATERIAL = 0,
    ARM_FACE_GROUP = 1,
    ARM_FACE_SMOOTHING_GROUPS = 2, /* bit mask per face, passed as int32 */
    ARM_FACE_MIRROR_SIDE = 3,      /* stage data only: 0 computed half, 1 mirrored half */
    ARM_FACE_PART = 4              /* stage data only: number of the part the face belongs to */
} arm_face_values;

/* Kinds of guide curves. Values are stable. */
typedef enum arm_curve_kind {
    ARM_CURVE_BAND = 0,      /* directions follow the curve in a band beside it */
    ARM_CURVE_EDGE_LOOP = 1, /* cut into the mesh: an edge loop follows the curve */
    ARM_CURVE_CREASE = 2,    /* like an edge loop; sharp in the result only for sharp lines of the preparation */
    ARM_CURVE_OFF = 3        /* ignored */
} arm_curve_kind;

/* Kinds of stage lines. Values are stable. */
typedef enum arm_line_kind {
    ARM_LINE_FEATURE = 0,    /* feature line of the preparation */
    ARM_LINE_GUIDE_CURVE = 1 /* guide curve of the host */
} arm_line_kind;

/* Settings block. The caller sets struct_size = sizeof(arm_settings); fields
 * beyond a smaller struct_size keep their defaults, so later versions can
 * append fields. Booleans are 0/1. */
typedef struct arm_settings {
    uint32_t struct_size;

    /* size and density */
    int32_t quad_count;          /* target number of quads (>= 1) */
    float quad_percentage;       /* > 0: quad count = input faces * percentage / 100 */
    float edge_length;           /* > 0: target edge length, takes precedence over the count */
    uint8_t exact_count;         /* keep the count strictly */
    uint8_t color_density;       /* density from the color attribute */
    uint8_t attribute_density;   /* density from the density attribute */
    float curvature;             /* 0..100: how strongly curvature densifies */
    float points_per_circle;     /* > 0 overrides the value derived from curvature */
    float density_ratio;         /* > 0 overrides the value derived from curvature */
    float max_quad_ratio;        /* largest side ratio */

    /* edge flow */
    uint8_t follow_materials;
    uint8_t follow_face_groups;
    uint8_t follow_smoothing_groups;
    uint8_t follow_sharp_edges;
    uint8_t follow_normal_splits;
    uint8_t detect_creases;
    uint8_t merge_cad_groups;
    uint8_t weld_seams;
    uint8_t weld_as_crease;
    uint8_t denoise;
    float crease_angle;          /* degrees */

    /* borders */
    uint8_t follow_borders;
    uint8_t freeze_borders;
    float border_density_factor;

    /* symmetry */
    uint8_t mirror_x;
    uint8_t mirror_y;
    uint8_t mirror_z;
    uint8_t mirror_at_center;    /* 1: plane through the bounding box center, 0: object origin */
    uint8_t detect_symmetry;
    uint8_t freeze_half;

    /* output */
    uint8_t split_at_creases;
    uint8_t creases_to_edge_crease;
    uint8_t clean_output;
    uint8_t keep_colors;

    /* computing */
    int32_t thread_count;        /* threads for the separate pieces of the mesh, 0 = all cores;
                                  * the result is the same for every count */
} arm_settings;

/* Progress callback: fraction 0..1 of the whole run and the current stage; always called on the
 * thread of arm_run. */
typedef void (*arm_progress_callback)(void* user, float fraction, int32_t stage);

/* Opaque handle owning all state of one remeshing session. */
typedef struct arm_context arm_context;

/* Version of the core library (major changes break the interface). */
ARM_API void arm_version(int32_t* major, int32_t* minor, int32_t* patch);

/* Version as "major.minor.patch"; static storage, never freed. */
ARM_API const char* arm_version_string(void);

/* Short English description of a status code; static storage. */
ARM_API const char* arm_status_text(arm_status status);

/* Creates a context. On success *out_context receives the handle. */
ARM_API arm_status arm_context_create(arm_context** out_context);

/* Destroys a context; null is allowed. */
ARM_API void arm_context_destroy(arm_context* context);

/* Message of the last failed call on this context, "" if none.
 * Valid until the next call on the same context. */
ARM_API const char* arm_context_last_error(const arm_context* context);

/* ---- remembered preparation ---- */

/* Opaque handle of a remembered preparation: the state after the stage ARM_STAGE_PREPARED of the
 * last run that computed it, with its key (a checksum of the input parts the stages up to it read
 * and every setting they read). A run on a context using the cache skips the stages
 * ARM_STAGE_INPUT and ARM_STAGE_PREPARED when the key matches, with the same result as without
 * the cache; otherwise it leaves its own preparation there. */
typedef struct arm_prep_cache arm_prep_cache;

ARM_API arm_status arm_prep_cache_create(arm_prep_cache** out_cache);
/* Gives up the host's hold; null is allowed. The memory is freed once no context uses it. */
ARM_API void arm_prep_cache_destroy(arm_prep_cache* cache);
/* Runs on this context use the cache from now on (null: none); the context holds it. */
ARM_API arm_status arm_context_use_prep_cache(arm_context* context, arm_prep_cache* cache);

/* ---- settings ---- */

/* Fills the block with the defaults; sets struct_size. */
ARM_API void arm_settings_default(arm_settings* settings);
ARM_API arm_status arm_settings_set(arm_context* context, const arm_settings* settings);

/* ---- input (object coordinates, corner order = face order) ---- */

/* Starts a new input and drops every earlier input part and stage result. */
ARM_API arm_status arm_input_mesh(arm_context* context, int32_t point_count, const float* points_xyz,
                                  int32_t face_count, const int32_t* face_sizes, int32_t corner_count,
                                  const int32_t* corner_points);
/* One value per face. */
ARM_API arm_status arm_input_face_values(arm_context* context, arm_face_values kind, const int32_t* values);
/* Sharp edges as point pairs; pairs without a mesh edge are ignored. */
ARM_API arm_status arm_input_sharp_edges(arm_context* context, int32_t edge_count, const int32_t* point_pairs);
/* One normal (xyz) per corner. share_equal != 0: equal normals share one index
 * (ascending order of the values); 0: every corner keeps its own normal. */
ARM_API arm_status arm_input_corner_normals(arm_context* context, const float* normals_xyz, int32_t share_equal);
/* One color (rgba, 0..1, sRGB) per corner; the last corner of a point decides its color. */
ARM_API arm_status arm_input_corner_colors(arm_context* context, const float* colors_rgba);
/* One density factor per point. */
ARM_API arm_status arm_input_point_density(arm_context* context, const float* values);
/* One UV pair per corner. */
ARM_API arm_status arm_input_uv_set(arm_context* context, const char* name, const float* uv);
/* One guide curve (after arm_input_mesh): its points (xyz, object coordinates of the
 * mesh) in curve order, its kind, on/off, width (< 0: automatic) and strength 0..100.
 * A closed curve repeats its first point at the end. */
ARM_API arm_status arm_input_curve(arm_context* context, int32_t point_count, const float* points_xyz,
                                   arm_curve_kind kind, int32_t active, float width, float strength);

/* ---- run ---- */

/* Runs until the given stage (inclusive). The callback may be null. ARM_STAGE_RESULT fails with
 * ARM_ERROR_NO_RESULT when no face is left. */
ARM_API arm_status arm_run(arm_context* context, arm_stage until, arm_progress_callback progress, void* user);
/* Asks a running arm_run to stop; safe from any thread. */
ARM_API void arm_cancel(arm_context* context);
/* Highest stage this build can compute (all lower ones are available too). */
ARM_API arm_stage arm_stage_last_available(void);

/* ---- data of the last computed stage ---- */

typedef struct arm_mesh_view {
    int32_t point_count;
    const float* points_xyz;
    int32_t face_count;
    const int32_t* face_starts;   /* face_count + 1 entries */
    int32_t corner_count;
    const int32_t* corner_points;
} arm_mesh_view;

ARM_API arm_status arm_stage_mesh(arm_context* context, arm_stage* out_stage, arm_mesh_view* out_mesh);
/* 1 when the run took the stages ARM_STAGE_INPUT and ARM_STAGE_PREPARED from the remembered
 * preparation, else 0. */
ARM_API arm_status arm_stage_preparation_reused(arm_context* context, int32_t* out_reused);
/* *out_values is null and *out_count 0 when the stage has no such values. */
ARM_API arm_status arm_stage_face_values(arm_context* context, arm_face_values kind, const int32_t** out_values,
                                         int32_t* out_count);
/* Sharp edges as point pairs (smaller point first), in edge order. In ARM_STAGE_RESULT every
 * other edge is smooth, and with creases_to_edge_crease these are the edges with full crease. */
ARM_API arm_status arm_stage_sharp_edges(arm_context* context, const int32_t** out_point_pairs,
                                         int32_t* out_edge_count);
/* Lines of the stage (from ARM_STAGE_FEATURE_LINES on): line i runs over the points
 * out_points[out_starts[i] .. out_starts[i + 1]) of the stage mesh; *out_line_count 0 when absent. */
ARM_API arm_status arm_stage_lines(arm_context* context, const int32_t** out_starts, const int32_t** out_points,
                                   int32_t* out_line_count);
/* Kind (arm_line_kind) of every stage line, in line order; *out_count 0 when absent. */
ARM_API arm_status arm_stage_line_kinds(arm_context* context, const int32_t** out_kinds, int32_t* out_count);
/* Parts of the stage (from ARM_STAGE_PARTS on): per part its quad budget and 1 when it is
 * skipped (too small), else 0, as pairs; *out_part_count 0 when absent. */
ARM_API arm_status arm_stage_parts(arm_context* context, const int32_t** out_budget_skipped,
                                   int32_t* out_part_count);
/* Normal index per corner and the normals; *out_normal_count 0 when absent. */
ARM_API arm_status arm_stage_normals(arm_context* context, const int32_t** out_corner_normal,
                                     const float** out_normals_xyz, int32_t* out_normal_count);
/* Color per point (rgba bytes); *out_count 0 when absent. */
ARM_API arm_status arm_stage_point_colors(arm_context* context, const uint8_t** out_rgba, int32_t* out_count);
ARM_API arm_status arm_stage_point_density(arm_context* context, const float** out_values, int32_t* out_count);
/* Curvature of the stage (ARM_STAGE_CURVATURE): per face the direction of least curvature
 * (xyz, along a limb), zero where none was found; *out_face_count 0 when absent. */
ARM_API arm_status arm_stage_curvature_directions(arm_context* context, const float** out_directions_xyz,
                                                  int32_t* out_face_count);
/* Cross field of the stage (ARM_STAGE_CROSS_FIELD): per face one arm of its cross (xyz, unit
 * length; the other arms follow by quarter turns about the face normal), zero where the face
 * has none; *out_face_count 0 when absent. */
ARM_API arm_status arm_stage_cross_directions(arm_context* context, const float** out_directions_xyz,
                                              int32_t* out_face_count);
/* Singular points of the cross field as pairs: the point of the stage mesh and the turn of the
 * crosses around it in quarter turns (1: three quads will meet there, -1: five, 2: half a turn;
 * 4 is added where the crosses make a whole turn); *out_count 0 when absent. */
ARM_API arm_status arm_stage_singularities(arm_context* context, const int32_t** out_point_index,
                                           int32_t* out_count);
/* Cuts of the stage (ARM_STAGE_CUTS and ARM_STAGE_PARAMETERIZATION): the edges along which
 * every part with a cross field is cut open into a disc, as point pairs of the stage mesh (a
 * part without cross field has no cuts); *out_edge_count 0 when absent. */
ARM_API arm_status arm_stage_cut_edges(arm_context* context, const int32_t** out_point_pairs,
                                       int32_t* out_edge_count);
/* Size field of the stage (ARM_STAGE_SIZE_FIELD): per face the side lengths the quads are to
 * get there, along the two directions of its cross (pairs; zero where the face has none);
 * *out_face_count 0 when absent. */
ARM_API arm_status arm_stage_size_field(arm_context* context, const float** out_sizes_uv,
                                        int32_t* out_face_count);
/* Per part of that stage four values: its target edge length, the smallest and the largest
 * size curvature may give a quad side, and the number of quads the size field amounts to;
 * *out_part_count 0 when absent. */
ARM_API arm_status arm_stage_size_limits(arm_context* context, const float** out_limits,
                                         int32_t* out_part_count);
/* Parameterization of the stage (ARM_STAGE_PARAMETERIZATION): per corner of the stage mesh its
 * two grid coordinates (pairs, in corner order); the quads of the result are the unit cells of
 * this grid. Corners at one point differ where the part is cut open. Zero for the corners of a
 * part without solution; *out_corner_count 0 when absent. */
ARM_API arm_status arm_stage_grid_coordinates(arm_context* context, const float** out_corner_uv,
                                              int32_t* out_corner_count);
/* Per part of that stage 1 when it has grid coordinates, else 0; *out_part_count 0 when absent. */
ARM_API arm_status arm_stage_grid_parts(arm_context* context, const int32_t** out_solved,
                                        int32_t* out_part_count);
ARM_API arm_status arm_stage_uv_set_count(arm_context* context, int32_t* out_count);
/* UV set `index`: name, UV point per corner, coordinates. */
ARM_API arm_status arm_stage_uv_set(arm_context* context, int32_t index, const char** out_name,
                                    const int32_t** out_corner_uv, const float** out_uv, int32_t* out_uv_count);

#ifdef __cplusplus
}
#endif

#endif /* AREKO_REMESHER_H */
