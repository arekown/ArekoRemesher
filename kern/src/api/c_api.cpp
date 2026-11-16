// Thin C shell around the core classes. No logic lives here: each function
// unwraps the handle, forwards to a class and lets guarded() map failures.

#include "areko_remesher.h"
#include "context.h"
#include "run/pipeline.h"
#include "run/prep_cache.h"
#include "settings/settings.h"
#include "status.h"
#include "version.h"

#include <memory>
#include <span>

using areko::remesher::api::Context;
using areko::remesher::api::guarded;
using areko::remesher::api::StageView;
using areko::remesher::api::StatusError;
using areko::remesher::api::Version;

struct arm_context {
    Context impl;
};

// The host's hold on a cache; every context using it holds it as well.
struct arm_prep_cache {
    std::shared_ptr<areko::remesher::run::PrepCache> impl;
};

namespace {

Context& unwrap(arm_context* context)
{
    if (context == nullptr) {
        throw StatusError(ARM_ERROR_INVALID_ARGUMENT, "context is null");
    }
    return context->impl;
}

Context* errorTarget(arm_context* context)
{
    return context != nullptr ? &context->impl : nullptr;
}

template <class T>
std::span<const T> items(const T* data, long long count, const char* what)
{
    if (count < 0) {
        throw StatusError(ARM_ERROR_INVALID_ARGUMENT, std::string(what) + ": negative count");
    }
    if (count > 0 && data == nullptr) {
        throw StatusError(ARM_ERROR_INVALID_ARGUMENT, std::string(what) + " is null");
    }
    return {data, static_cast<std::size_t>(count)};
}

template <class T>
void require(T* out, const char* what)
{
    if (out == nullptr) {
        throw StatusError(ARM_ERROR_INVALID_ARGUMENT, std::string(what) + " is null");
    }
}

} // namespace

extern "C" {

ARM_API void arm_version(int32_t* major, int32_t* minor, int32_t* patch)
{
    if (major != nullptr) {
        *major = Version::major;
    }
    if (minor != nullptr) {
        *minor = Version::minor;
    }
    if (patch != nullptr) {
        *patch = Version::patch;
    }
}

ARM_API const char* arm_version_string(void)
{
    return Version::text;
}

ARM_API const char* arm_status_text(arm_status status)
{
    return areko::remesher::api::statusText(status);
}

ARM_API arm_status arm_context_create(arm_context** out_context)
{
    return guarded(nullptr, [out_context] {
        if (out_context == nullptr) {
            throw StatusError(ARM_ERROR_INVALID_ARGUMENT, "out_context is null");
        }
        *out_context = nullptr;
        *out_context = std::make_unique<arm_context>().release();
    });
}

ARM_API void arm_context_destroy(arm_context* context)
{
    std::unique_ptr<arm_context> owner(context);
}

ARM_API const char* arm_context_last_error(const arm_context* context)
{
    return context != nullptr ? context->impl.lastError().c_str() : "";
}

ARM_API arm_status arm_prep_cache_create(arm_prep_cache** out_cache)
{
    return guarded(nullptr, [out_cache] {
        if (out_cache == nullptr) {
            throw StatusError(ARM_ERROR_INVALID_ARGUMENT, "out_cache is null");
        }
        *out_cache = nullptr;
        auto cache = std::make_unique<arm_prep_cache>();
        cache->impl = std::make_shared<areko::remesher::run::PrepCache>();
        *out_cache = cache.release();
    });
}

ARM_API void arm_prep_cache_destroy(arm_prep_cache* cache)
{
    std::unique_ptr<arm_prep_cache> owner(cache);
}

ARM_API arm_status arm_context_use_prep_cache(arm_context* context, arm_prep_cache* cache)
{
    return guarded(errorTarget(context), [&] {
        unwrap(context).usePrepCache(cache != nullptr ? cache->impl : nullptr);
    });
}

ARM_API void arm_settings_default(arm_settings* settings)
{
    if (settings != nullptr) {
        areko::remesher::settings::Settings::fillDefaults(*settings);
    }
}

ARM_API arm_status arm_settings_set(arm_context* context, const arm_settings* settings)
{
    return guarded(errorTarget(context), [&] {
        require(settings, "settings");
        unwrap(context).settings().assign(*settings);
    });
}

ARM_API arm_status arm_input_mesh(arm_context* context, int32_t point_count, const float* points_xyz,
                                  int32_t face_count, const int32_t* face_sizes, int32_t corner_count,
                                  const int32_t* corner_points)
{
    return guarded(errorTarget(context), [&] {
        auto& input = unwrap(context).startInput();
        input.setMesh(items(points_xyz, 3LL * point_count, "points"), items(face_sizes, face_count, "face sizes"),
                      items(corner_points, corner_count, "corner points"));
    });
}

ARM_API arm_status arm_input_face_values(arm_context* context, arm_face_values kind, const int32_t* values)
{
    return guarded(errorTarget(context), [&] {
        auto& input = unwrap(context).editInput();
        input.setFaceValues(kind, items(values, input.faceCount(), "face values"));
    });
}

ARM_API arm_status arm_input_sharp_edges(arm_context* context, int32_t edge_count, const int32_t* point_pairs)
{
    return guarded(errorTarget(context), [&] {
        unwrap(context).editInput().setSharpEdges(items(point_pairs, 2LL * edge_count, "point pairs"));
    });
}

ARM_API arm_status arm_input_corner_normals(arm_context* context, const float* normals_xyz, int32_t share_equal)
{
    return guarded(errorTarget(context), [&] {
        auto& input = unwrap(context).editInput();
        input.setCornerNormals(items(normals_xyz, 3LL * input.cornerCount(), "normals"), share_equal != 0);
    });
}

ARM_API arm_status arm_input_corner_colors(arm_context* context, const float* colors_rgba)
{
    return guarded(errorTarget(context), [&] {
        auto& input = unwrap(context).editInput();
        input.setCornerColors(items(colors_rgba, 4LL * input.cornerCount(), "colors"));
    });
}

ARM_API arm_status arm_input_point_density(arm_context* context, const float* values)
{
    return guarded(errorTarget(context), [&] {
        auto& input = unwrap(context).editInput();
        input.setPointDensity(items(values, input.pointCount(), "density"));
    });
}

ARM_API arm_status arm_input_uv_set(arm_context* context, const char* name, const float* uv)
{
    return guarded(errorTarget(context), [&] {
        auto& input = unwrap(context).editInput();
        input.addUvSet(name != nullptr ? name : "", items(uv, 2LL * input.cornerCount(), "uv"));
    });
}

ARM_API arm_status arm_input_curve(arm_context* context, int32_t point_count, const float* points_xyz,
                                   arm_curve_kind kind, int32_t active, float width, float strength)
{
    return guarded(errorTarget(context), [&] {
        unwrap(context).editInput().addCurve(items(points_xyz, 3LL * point_count, "points_xyz"),
                                             static_cast<int>(kind), active != 0, width, strength);
    });
}

ARM_API arm_status arm_run(arm_context* context, arm_stage until, arm_progress_callback progress, void* user)
{
    return guarded(errorTarget(context), [&] { unwrap(context).run(until, progress, user); });
}

ARM_API void arm_cancel(arm_context* context)
{
    if (context != nullptr) {
        context->impl.requestCancel();
    }
}

ARM_API arm_stage arm_stage_last_available(void)
{
    return areko::remesher::run::Pipeline::lastAvailable;
}

ARM_API arm_status arm_stage_mesh(arm_context* context, arm_stage* out_stage, arm_mesh_view* out_mesh)
{
    return guarded(errorTarget(context), [&] {
        require(out_mesh, "out_mesh");
        const StageView& stage = unwrap(context).stage();
        if (out_stage != nullptr) {
            *out_stage = stage.stage();
        }
        *out_mesh = stage.meshView();
    });
}

ARM_API arm_status arm_stage_preparation_reused(arm_context* context, int32_t* out_reused)
{
    return guarded(errorTarget(context), [&] {
        require(out_reused, "out_reused");
        *out_reused = unwrap(context).preparationReused() ? 1 : 0;
    });
}

ARM_API arm_status arm_stage_face_values(arm_context* context, arm_face_values kind, const int32_t** out_values,
                                         int32_t* out_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_values, "out_values");
        require(out_count, "out_count");
        const auto& values = unwrap(context).stage().faceValues(kind);
        *out_values = values ? values->data() : nullptr;
        *out_count = values ? static_cast<int32_t>(values->size()) : 0;
    });
}

ARM_API arm_status arm_stage_sharp_edges(arm_context* context, const int32_t** out_point_pairs,
                                         int32_t* out_edge_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_point_pairs, "out_point_pairs");
        require(out_edge_count, "out_edge_count");
        const auto& pairs = unwrap(context).stage().sharpEdgePairs();
        *out_point_pairs = pairs.empty() ? nullptr : pairs.data();
        *out_edge_count = static_cast<int32_t>(pairs.size() / 2);
    });
}

ARM_API arm_status arm_stage_lines(arm_context* context, const int32_t** out_starts, const int32_t** out_points,
                                   int32_t* out_line_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_starts, "out_starts");
        require(out_points, "out_points");
        require(out_line_count, "out_line_count");
        const StageView& stage = unwrap(context).stage();
        const auto& starts = stage.lineStarts();
        const auto& points = stage.linePoints();
        const auto count = static_cast<int32_t>(starts.size()) - 1;
        *out_starts = count > 0 ? starts.data() : nullptr;
        *out_points = points.empty() ? nullptr : points.data();
        *out_line_count = count;
    });
}

ARM_API arm_status arm_stage_line_kinds(arm_context* context, const int32_t** out_kinds, int32_t* out_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_kinds, "out_kinds");
        require(out_count, "out_count");
        const auto& kinds = unwrap(context).stage().lineKinds();
        *out_kinds = kinds.empty() ? nullptr : kinds.data();
        *out_count = static_cast<int32_t>(kinds.size());
    });
}

ARM_API arm_status arm_stage_parts(arm_context* context, const int32_t** out_budget_skipped, int32_t* out_part_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_budget_skipped, "out_budget_skipped");
        require(out_part_count, "out_part_count");
        const auto& values = unwrap(context).stage().partValues();
        *out_budget_skipped = values.empty() ? nullptr : values.data();
        *out_part_count = static_cast<int32_t>(values.size() / 2);
    });
}

ARM_API arm_status arm_stage_normals(arm_context* context, const int32_t** out_corner_normal,
                                     const float** out_normals_xyz, int32_t* out_normal_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_corner_normal, "out_corner_normal");
        require(out_normals_xyz, "out_normals_xyz");
        require(out_normal_count, "out_normal_count");
        const StageView& stage = unwrap(context).stage();
        const bool present = !stage.cornerNormal().empty();
        *out_corner_normal = present ? stage.cornerNormal().data() : nullptr;
        *out_normals_xyz = present ? stage.normals().data() : nullptr;
        *out_normal_count = present ? static_cast<int32_t>(stage.normals().size() / 3) : 0;
    });
}

ARM_API arm_status arm_stage_point_colors(arm_context* context, const uint8_t** out_rgba, int32_t* out_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_rgba, "out_rgba");
        require(out_count, "out_count");
        const auto& colors = unwrap(context).stage().colors();
        *out_rgba = colors.empty() ? nullptr : colors.data();
        *out_count = static_cast<int32_t>(colors.size() / 4);
    });
}

ARM_API arm_status arm_stage_point_density(arm_context* context, const float** out_values, int32_t* out_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_values, "out_values");
        require(out_count, "out_count");
        const auto& density = unwrap(context).stage().density();
        *out_values = density.empty() ? nullptr : density.data();
        *out_count = static_cast<int32_t>(density.size());
    });
}

ARM_API arm_status arm_stage_curvature_directions(arm_context* context, const float** out_directions_xyz,
                                                  int32_t* out_face_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_directions_xyz, "out_directions_xyz");
        require(out_face_count, "out_face_count");
        const auto& directions = unwrap(context).stage().faceDirections();
        *out_directions_xyz = directions.empty() ? nullptr : directions.data();
        *out_face_count = static_cast<int32_t>(directions.size() / 3);
    });
}

ARM_API arm_status arm_stage_cross_directions(arm_context* context, const float** out_directions_xyz,
                                              int32_t* out_face_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_directions_xyz, "out_directions_xyz");
        require(out_face_count, "out_face_count");
        const auto& directions = unwrap(context).stage().crossDirections();
        *out_directions_xyz = directions.empty() ? nullptr : directions.data();
        *out_face_count = static_cast<int32_t>(directions.size() / 3);
    });
}

ARM_API arm_status arm_stage_singularities(arm_context* context, const int32_t** out_point_index, int32_t* out_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_point_index, "out_point_index");
        require(out_count, "out_count");
        const auto& pairs = unwrap(context).stage().singularities();
        *out_point_index = pairs.empty() ? nullptr : pairs.data();
        *out_count = static_cast<int32_t>(pairs.size() / 2);
    });
}

ARM_API arm_status arm_stage_cut_edges(arm_context* context, const int32_t** out_point_pairs, int32_t* out_edge_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_point_pairs, "out_point_pairs");
        require(out_edge_count, "out_edge_count");
        const auto& pairs = unwrap(context).stage().cutEdgePairs();
        *out_point_pairs = pairs.empty() ? nullptr : pairs.data();
        *out_edge_count = static_cast<int32_t>(pairs.size() / 2);
    });
}

ARM_API arm_status arm_stage_size_field(arm_context* context, const float** out_sizes_uv, int32_t* out_face_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_sizes_uv, "out_sizes_uv");
        require(out_face_count, "out_face_count");
        const auto& sizes = unwrap(context).stage().faceSizes();
        *out_sizes_uv = sizes.empty() ? nullptr : sizes.data();
        *out_face_count = static_cast<int32_t>(sizes.size() / 2);
    });
}

ARM_API arm_status arm_stage_size_limits(arm_context* context, const float** out_limits, int32_t* out_part_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_limits, "out_limits");
        require(out_part_count, "out_part_count");
        const auto& limits = unwrap(context).stage().partSizes();
        *out_limits = limits.empty() ? nullptr : limits.data();
        *out_part_count = static_cast<int32_t>(limits.size() / 4);
    });
}

ARM_API arm_status arm_stage_grid_coordinates(arm_context* context, const float** out_corner_uv,
                                              int32_t* out_corner_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_corner_uv, "out_corner_uv");
        require(out_corner_count, "out_corner_count");
        const auto& coordinates = unwrap(context).stage().gridCoordinates();
        *out_corner_uv = coordinates.empty() ? nullptr : coordinates.data();
        *out_corner_count = static_cast<int32_t>(coordinates.size() / 2);
    });
}

ARM_API arm_status arm_stage_grid_parts(arm_context* context, const int32_t** out_solved, int32_t* out_part_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_solved, "out_solved");
        require(out_part_count, "out_part_count");
        const auto& parts = unwrap(context).stage().gridParts();
        *out_solved = parts.empty() ? nullptr : parts.data();
        *out_part_count = static_cast<int32_t>(parts.size());
    });
}

ARM_API arm_status arm_stage_uv_set_count(arm_context* context, int32_t* out_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_count, "out_count");
        *out_count = static_cast<int32_t>(unwrap(context).stage().uvSets().size());
    });
}

ARM_API arm_status arm_stage_uv_set(arm_context* context, int32_t index, const char** out_name,
                                    const int32_t** out_corner_uv, const float** out_uv, int32_t* out_uv_count)
{
    return guarded(errorTarget(context), [&] {
        require(out_name, "out_name");
        require(out_corner_uv, "out_corner_uv");
        require(out_uv, "out_uv");
        require(out_uv_count, "out_uv_count");
        const auto& sets = unwrap(context).stage().uvSets();
        if (index < 0 || static_cast<std::size_t>(index) >= sets.size()) {
            throw StatusError(ARM_ERROR_INVALID_ARGUMENT, "no such UV set");
        }
        const auto& set = sets[static_cast<std::size_t>(index)];
        *out_name = set.name.c_str();
        *out_corner_uv = set.cornerUv.data();
        *out_uv = set.coords.data();
        *out_uv_count = static_cast<int32_t>(set.coords.size() / 2);
    });
}

} // extern "C"
