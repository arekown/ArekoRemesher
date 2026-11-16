# Binding of the core library (libarekoremesher.so / arekoremesher.dll / libarekoremesher.dylib) via ctypes.
# All C signatures are declared in one place (CoreLibrary._declare).

from __future__ import annotations

import ctypes
import sys
import threading
import traceback
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ARM_OK = 0
ARM_CANCELLED = 6

_i32 = ctypes.c_int32
_i32p = ctypes.POINTER(ctypes.c_int32)
_f32p = ctypes.POINTER(ctypes.c_float)
_u8p = ctypes.POINTER(ctypes.c_uint8)
_u8 = ctypes.c_uint8


class CoreSettings(ctypes.Structure):
    """Mirror of arm_settings (same field order as the C header)."""

    _fields_ = [
        ("struct_size", ctypes.c_uint32),
        ("quad_count", _i32),
        ("quad_percentage", ctypes.c_float),
        ("edge_length", ctypes.c_float),
        ("exact_count", _u8),
        ("color_density", _u8),
        ("attribute_density", _u8),
        ("curvature", ctypes.c_float),
        ("points_per_circle", ctypes.c_float),
        ("density_ratio", ctypes.c_float),
        ("max_quad_ratio", ctypes.c_float),
        ("follow_materials", _u8),
        ("follow_face_groups", _u8),
        ("follow_smoothing_groups", _u8),
        ("follow_sharp_edges", _u8),
        ("follow_normal_splits", _u8),
        ("detect_creases", _u8),
        ("merge_cad_groups", _u8),
        ("weld_seams", _u8),
        ("weld_as_crease", _u8),
        ("denoise", _u8),
        ("crease_angle", ctypes.c_float),
        ("follow_borders", _u8),
        ("freeze_borders", _u8),
        ("border_density_factor", ctypes.c_float),
        ("mirror_x", _u8),
        ("mirror_y", _u8),
        ("mirror_z", _u8),
        ("mirror_at_center", _u8),
        ("detect_symmetry", _u8),
        ("freeze_half", _u8),
        ("split_at_creases", _u8),
        ("creases_to_edge_crease", _u8),
        ("clean_output", _u8),
        ("keep_colors", _u8),
        ("thread_count", _i32),
    ]

    @property
    def field_names(self) -> tuple[str, ...]:
        return tuple(name for name, _ in self._fields_ if name != "struct_size")


class _MeshView(ctypes.Structure):
    _fields_ = [
        ("point_count", _i32),
        ("points_xyz", _f32p),
        ("face_count", _i32),
        ("face_starts", _i32p),
        ("corner_count", _i32),
        ("corner_points", _i32p),
    ]


_PROGRESS = ctypes.CFUNCTYPE(None, ctypes.c_void_p, ctypes.c_float, _i32)

# Face value kinds (arm_face_values).
FACE_MATERIAL = 0
FACE_GROUP = 1
FACE_SMOOTHING_GROUPS = 2
FACE_MIRROR_SIDE = 3  # stage data only: 0 computed half, 1 mirrored half
FACE_PART = 4  # stage data only: part number

# Guide curve kinds (arm_curve_kind) and stage line kinds (arm_line_kind).
CURVE_BAND, CURVE_EDGE_LOOP, CURVE_CREASE, CURVE_OFF = 0, 1, 2, 3
LINE_FEATURE, LINE_GUIDE_CURVE = 0, 1

# Status codes (arm_status) -> message shown to the user.
STATUS_MESSAGES = {
    1: "Invalid input",
    2: "Out of memory",
    3: "Internal error",
    4: "No mesh",
    5: "No result",
    6: "Cancelled",
    7: "Invalid settings",
    8: "Stage not available yet",
}


class CoreError(RuntimeError):
    """Raised when the core library cannot be used."""


class CoreRunError(CoreError):
    """A core call failed with a status code."""

    def __init__(self, status: int, detail: str):
        super().__init__(f"{STATUS_MESSAGES.get(status, 'Internal error')}: {detail}")
        self.status = status
        self.detail = detail

    @property
    def message(self) -> str:
        return STATUS_MESSAGES.get(self.status, "Internal error")


def _ptr(array: np.ndarray, ctype):
    return array.ctypes.data_as(ctypes.POINTER(ctype))


def _copy(pointer, count: int, dtype, width: int = 1) -> np.ndarray:
    if count <= 0 or not pointer:
        return np.empty((0, width) if width > 1 else 0, dtype=dtype)
    flat = np.ctypeslib.as_array(pointer, shape=(count * width,)).astype(dtype, copy=True)
    return flat.reshape(-1, width) if width > 1 else flat


@dataclass
class CurveInput:
    """One guide curve: points (object coordinates of the mesh) in curve order."""

    points: np.ndarray          # float32 (n, 3)
    kind: int = CURVE_EDGE_LOOP
    active: bool = True
    width: float = -1.0         # < 0: automatic
    strength: float = 50.0      # 0..100


@dataclass
class MeshInput:
    """Raw input of one run (object coordinates, corners in face order)."""

    points: np.ndarray                      # float32 (n, 3)
    face_sizes: np.ndarray                  # int32 (f,)
    corner_points: np.ndarray               # int32 (c,)
    materials: np.ndarray | None = None     # int32 (f,)
    face_groups: np.ndarray | None = None   # int32 (f,)
    smoothing_groups: np.ndarray | None = None  # int32 bit masks (f,)
    sharp_edges: np.ndarray | None = None   # int32 (k, 2) point pairs
    corner_normals: np.ndarray | None = None    # float32 (c, 3)
    share_normals: bool = True                  # equal corner normals share one index
    corner_colors: np.ndarray | None = None     # float32 (c, 4), sRGB 0..1
    point_density: np.ndarray | None = None     # float32 (n,)
    uv_sets: list[tuple[str, np.ndarray]] = field(default_factory=list)  # float32 (c, 2)
    curves: list[CurveInput] = field(default_factory=list)


@dataclass
class StageData:
    """Copy of the data of one computed stage."""

    stage: int
    points: np.ndarray
    face_starts: np.ndarray
    corner_points: np.ndarray
    materials: np.ndarray | None
    face_groups: np.ndarray | None
    smoothing_groups: np.ndarray | None
    mirror_sides: np.ndarray | None
    face_parts: np.ndarray | None
    part_values: np.ndarray        # int32 (parts, 2): quad budget, 1 = skipped
    sharp_edges: np.ndarray
    line_starts: np.ndarray        # int32 (lines + 1,): line i runs over line_points[starts[i]:starts[i + 1]]
    line_points: np.ndarray        # int32 point numbers of the stage mesh
    line_kinds: np.ndarray         # int32 per line (LINE_*); empty when the stage has none
    corner_normal: np.ndarray | None
    normals: np.ndarray | None
    colors: np.ndarray | None
    density: np.ndarray | None
    uv_sets: list[tuple[str, np.ndarray, np.ndarray]]
    curvature_directions: np.ndarray | None = None  # float32 (faces, 3): direction of least curvature, zero = none
    cross_directions: np.ndarray | None = None  # float32 (faces, 3): one arm of the cross, zero = none
    singularities: np.ndarray | None = None  # int32 (n, 2): point of the stage mesh, turn of the crosses in quarter turns
    cut_edges: np.ndarray | None = None  # int32 (n, 2): point pairs of the stage mesh along which the parts are cut open
    size_field: np.ndarray | None = None  # float32 (faces, 2): side lengths of the quads along the two cross directions
    size_limits: np.ndarray | None = None  # float32 (parts, 4): edge length, smallest size, largest size, estimated quads
    grid_uv: np.ndarray | None = None  # float32 (corners, 2): grid coordinates; the quads are the unit cells
    grid_parts: np.ndarray | None = None  # int32 (parts,): 1 = the part has grid coordinates


class CoreLibrary:
    """Loaded core library with typed entry points."""

    if sys.platform == "win32":
        FILE_NAME = "arekoremesher.dll"
    elif sys.platform == "darwin":
        FILE_NAME = "libarekoremesher.dylib"
    else:
        FILE_NAME = "libarekoremesher.so"
    SUPPORTED_MAJOR = 0

    def __init__(self, path: Path):
        try:
            self._lib = ctypes.CDLL(str(path))
        except OSError as error:
            raise CoreError(f"cannot load {path.name}: {error}") from error
        self._declare()
        major = self.version_tuple[0]
        if major != self.SUPPORTED_MAJOR:
            raise CoreError(f"core version {self.version} is not supported")
        if self.default_settings().struct_size != ctypes.sizeof(CoreSettings):
            raise CoreError("settings block of the core does not match this add-on")

    @classmethod
    def from_addon(cls) -> CoreLibrary:
        return cls(Path(__file__).resolve().parent / "lib" / cls.FILE_NAME)

    def _declare(self) -> None:
        lib = self._lib
        vp = ctypes.c_void_p

        def sig(name, args, res=ctypes.c_int):
            fn = getattr(lib, name)
            fn.argtypes = args
            fn.restype = res

        sig("arm_version", [_i32p, _i32p, _i32p], None)
        sig("arm_version_string", [], ctypes.c_char_p)
        sig("arm_status_text", [ctypes.c_int], ctypes.c_char_p)
        sig("arm_context_create", [ctypes.POINTER(vp)])
        sig("arm_context_destroy", [vp], None)
        sig("arm_context_last_error", [vp], ctypes.c_char_p)
        sig("arm_prep_cache_create", [ctypes.POINTER(vp)])
        sig("arm_prep_cache_destroy", [vp], None)
        sig("arm_context_use_prep_cache", [vp, vp])
        sig("arm_settings_default", [ctypes.POINTER(CoreSettings)], None)
        sig("arm_settings_set", [vp, ctypes.POINTER(CoreSettings)])
        sig("arm_input_mesh", [vp, _i32, _f32p, _i32, _i32p, _i32, _i32p])
        sig("arm_input_face_values", [vp, ctypes.c_int, _i32p])
        sig("arm_input_sharp_edges", [vp, _i32, _i32p])
        sig("arm_input_corner_normals", [vp, _f32p, _i32])
        sig("arm_input_corner_colors", [vp, _f32p])
        sig("arm_input_point_density", [vp, _f32p])
        sig("arm_input_uv_set", [vp, ctypes.c_char_p, _f32p])
        sig("arm_input_curve", [vp, _i32, _f32p, ctypes.c_int, _i32, ctypes.c_float, ctypes.c_float])
        sig("arm_run", [vp, ctypes.c_int, _PROGRESS, vp])
        sig("arm_cancel", [vp], None)
        sig("arm_stage_last_available", [], ctypes.c_int)
        sig("arm_stage_mesh", [vp, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(_MeshView)])
        sig("arm_stage_preparation_reused", [vp, _i32p])
        sig("arm_stage_face_values", [vp, ctypes.c_int, ctypes.POINTER(_i32p), _i32p])
        sig("arm_stage_sharp_edges", [vp, ctypes.POINTER(_i32p), _i32p])
        sig("arm_stage_lines", [vp, ctypes.POINTER(_i32p), ctypes.POINTER(_i32p), _i32p])
        sig("arm_stage_line_kinds", [vp, ctypes.POINTER(_i32p), _i32p])
        sig("arm_stage_parts", [vp, ctypes.POINTER(_i32p), _i32p])
        sig("arm_stage_normals", [vp, ctypes.POINTER(_i32p), ctypes.POINTER(_f32p), _i32p])
        sig("arm_stage_point_colors", [vp, ctypes.POINTER(_u8p), _i32p])
        sig("arm_stage_point_density", [vp, ctypes.POINTER(_f32p), _i32p])
        sig("arm_stage_curvature_directions", [vp, ctypes.POINTER(_f32p), _i32p])
        sig("arm_stage_cross_directions", [vp, ctypes.POINTER(_f32p), _i32p])
        sig("arm_stage_singularities", [vp, ctypes.POINTER(_i32p), _i32p])
        sig("arm_stage_cut_edges", [vp, ctypes.POINTER(_i32p), _i32p])
        sig("arm_stage_size_field", [vp, ctypes.POINTER(_f32p), _i32p])
        sig("arm_stage_size_limits", [vp, ctypes.POINTER(_f32p), _i32p])
        sig("arm_stage_grid_coordinates", [vp, ctypes.POINTER(_f32p), _i32p])
        sig("arm_stage_grid_parts", [vp, ctypes.POINTER(_i32p), _i32p])
        sig("arm_stage_uv_set_count", [vp, _i32p])
        sig("arm_stage_uv_set", [vp, _i32, ctypes.POINTER(ctypes.c_char_p), ctypes.POINTER(_i32p),
                                 ctypes.POINTER(_f32p), _i32p])

    @property
    def version_tuple(self) -> tuple[int, int, int]:
        parts = [ctypes.c_int32() for _ in range(3)]
        self._lib.arm_version(*(ctypes.byref(p) for p in parts))
        return tuple(p.value for p in parts)

    @property
    def version(self) -> str:
        return self._lib.arm_version_string().decode("utf-8")

    @property
    def last_stage(self) -> int:
        return int(self._lib.arm_stage_last_available())

    def status_text(self, status: int) -> str:
        return self._lib.arm_status_text(status).decode("utf-8")

    def default_settings(self) -> CoreSettings:
        settings = CoreSettings()
        self._lib.arm_settings_default(ctypes.byref(settings))
        return settings

    def create_context(self) -> CoreContext:
        handle = ctypes.c_void_p()
        status = self._lib.arm_context_create(ctypes.byref(handle))
        if status != ARM_OK or not handle.value:
            raise CoreError(f"context: {self.status_text(status)}")
        return CoreContext(self._lib, handle)

    def create_prep_cache(self) -> CorePrepCache:
        handle = ctypes.c_void_p()
        status = self._lib.arm_prep_cache_create(ctypes.byref(handle))
        if status != ARM_OK or not handle.value:
            raise CoreError(f"prep cache: {self.status_text(status)}")
        return CorePrepCache(self._lib, handle)


class CorePrepCache:
    """Holds one remembered preparation of the core; release with close().

    Contexts that use it hold it as well, so close() may come while a run still
    uses it: the core frees it after that run.
    """

    def __init__(self, lib: ctypes.CDLL, handle: ctypes.c_void_p):
        self._lib = lib
        self._handle = handle
        self._lock = threading.Lock()  # guards the handle between close() and held()

    def __del__(self):
        self.close()

    def close(self) -> None:
        with self._lock:
            if self._handle is not None and self._handle.value:
                self._lib.arm_prep_cache_destroy(self._handle)
            self._handle = None

    @contextmanager
    def held(self):
        """The handle (None once closed), kept open for the with-block."""
        with self._lock:
            yield self._handle


class CoreContext:
    """Owns one core context; release with close() or a with-block."""

    def __init__(self, lib: ctypes.CDLL, handle: ctypes.c_void_p):
        self._lib = lib
        self._handle = handle

    def __enter__(self) -> CoreContext:
        return self

    def __exit__(self, *_exc) -> None:
        self.close()

    def __del__(self):
        self.close()

    def close(self) -> None:
        if self._handle is not None and self._handle.value:
            self._lib.arm_context_destroy(self._handle)
        self._handle = None

    @property
    def last_error(self) -> str:
        if self._handle is None:
            return ""
        return self._lib.arm_context_last_error(self._handle).decode("utf-8")

    def _check(self, status: int) -> None:
        if status != ARM_OK:
            raise CoreRunError(status, self.last_error)

    def set_settings(self, settings: CoreSettings) -> None:
        self._check(self._lib.arm_settings_set(self._handle, ctypes.byref(settings)))

    def load_input(self, mesh: MeshInput) -> None:
        h, lib = self._handle, self._lib
        points = np.ascontiguousarray(mesh.points, dtype=np.float32).reshape(-1)
        sizes = np.ascontiguousarray(mesh.face_sizes, dtype=np.int32)
        corners = np.ascontiguousarray(mesh.corner_points, dtype=np.int32)
        self._check(lib.arm_input_mesh(h, len(points) // 3, _ptr(points, ctypes.c_float), len(sizes),
                                       _ptr(sizes, ctypes.c_int32), len(corners), _ptr(corners, ctypes.c_int32)))
        for kind, values in ((FACE_GROUP, mesh.face_groups), (FACE_MATERIAL, mesh.materials),
                             (FACE_SMOOTHING_GROUPS, mesh.smoothing_groups)):
            if values is not None:
                ints = np.ascontiguousarray(values, dtype=np.int32)
                self._check(lib.arm_input_face_values(h, kind, _ptr(ints, ctypes.c_int32)))
        if mesh.sharp_edges is not None and len(mesh.sharp_edges):
            pairs = np.ascontiguousarray(mesh.sharp_edges, dtype=np.int32).reshape(-1)
            self._check(lib.arm_input_sharp_edges(h, len(pairs) // 2, _ptr(pairs, ctypes.c_int32)))
        for fn, values, extra in ((lib.arm_input_corner_normals, mesh.corner_normals, (int(mesh.share_normals),)),
                                  (lib.arm_input_corner_colors, mesh.corner_colors, ()),
                                  (lib.arm_input_point_density, mesh.point_density, ())):
            if values is not None:
                floats = np.ascontiguousarray(values, dtype=np.float32).reshape(-1)
                self._check(fn(h, _ptr(floats, ctypes.c_float), *extra))
        for name, uv in mesh.uv_sets:
            floats = np.ascontiguousarray(uv, dtype=np.float32).reshape(-1)
            self._check(lib.arm_input_uv_set(h, name.encode("utf-8"), _ptr(floats, ctypes.c_float)))
        for curve in mesh.curves:
            floats = np.ascontiguousarray(curve.points, dtype=np.float32).reshape(-1)
            self._check(lib.arm_input_curve(h, len(floats) // 3, _ptr(floats, ctypes.c_float), int(curve.kind),
                                            int(bool(curve.active)), float(curve.width), float(curve.strength)))

    def run(self, stage: int, progress=None) -> None:
        """Runs until `stage`; progress(fraction, stage) is called from the thread of the run."""
        callback = _PROGRESS(lambda _user, fraction, current: progress(fraction, current)) if progress \
            else _PROGRESS()
        self._check(self._lib.arm_run(self._handle, stage, callback, None))

    def cancel(self) -> None:
        if self._handle is not None:
            self._lib.arm_cancel(self._handle)

    def use_prep_cache(self, cache: CorePrepCache) -> None:
        """Runs on this context use `cache` from now on (none once it is closed)."""
        with cache.held() as handle:
            if handle is not None:
                self._check(self._lib.arm_context_use_prep_cache(self._handle, handle))

    @property
    def preparation_reused(self) -> bool:
        """The last run took its preparation from the remembered one."""
        reused = _i32()
        self._check(self._lib.arm_stage_preparation_reused(self._handle, ctypes.byref(reused)))
        return reused.value == 1

    def stage_data(self) -> StageData:
        h, lib = self._handle, self._lib
        view = _MeshView()
        stage = ctypes.c_int()
        self._check(lib.arm_stage_mesh(h, ctypes.byref(stage), ctypes.byref(view)))

        def face_values(kind):
            ptr, count = _i32p(), _i32()
            self._check(lib.arm_stage_face_values(h, kind, ctypes.byref(ptr), ctypes.byref(count)))
            return _copy(ptr, count.value, np.int32) if count.value else None

        pairs, edge_count = _i32p(), _i32()
        self._check(lib.arm_stage_sharp_edges(h, ctypes.byref(pairs), ctypes.byref(edge_count)))
        line_starts, line_points, line_count = _i32p(), _i32p(), _i32()
        self._check(lib.arm_stage_lines(h, ctypes.byref(line_starts), ctypes.byref(line_points),
                                        ctypes.byref(line_count)))
        starts = _copy(line_starts, line_count.value + 1, np.int32) if line_count.value > 0 \
            else np.zeros(1, dtype=np.int32)
        kinds, kind_count = _i32p(), _i32()
        self._check(lib.arm_stage_line_kinds(h, ctypes.byref(kinds), ctypes.byref(kind_count)))
        parts, part_count = _i32p(), _i32()
        self._check(lib.arm_stage_parts(h, ctypes.byref(parts), ctypes.byref(part_count)))
        corner_normal, normals, normal_count = _i32p(), _f32p(), _i32()
        self._check(lib.arm_stage_normals(h, ctypes.byref(corner_normal), ctypes.byref(normals),
                                          ctypes.byref(normal_count)))
        colors, color_count = _u8p(), _i32()
        self._check(lib.arm_stage_point_colors(h, ctypes.byref(colors), ctypes.byref(color_count)))
        density, density_count = _f32p(), _i32()
        self._check(lib.arm_stage_point_density(h, ctypes.byref(density), ctypes.byref(density_count)))
        directions, direction_count = _f32p(), _i32()
        self._check(lib.arm_stage_curvature_directions(h, ctypes.byref(directions), ctypes.byref(direction_count)))
        crosses, cross_count = _f32p(), _i32()
        self._check(lib.arm_stage_cross_directions(h, ctypes.byref(crosses), ctypes.byref(cross_count)))
        singular, singular_count = _i32p(), _i32()
        self._check(lib.arm_stage_singularities(h, ctypes.byref(singular), ctypes.byref(singular_count)))
        cuts, cut_count = _i32p(), _i32()
        self._check(lib.arm_stage_cut_edges(h, ctypes.byref(cuts), ctypes.byref(cut_count)))
        sizes, size_count = _f32p(), _i32()
        self._check(lib.arm_stage_size_field(h, ctypes.byref(sizes), ctypes.byref(size_count)))
        limits, limit_count = _f32p(), _i32()
        self._check(lib.arm_stage_size_limits(h, ctypes.byref(limits), ctypes.byref(limit_count)))
        grid, grid_count = _f32p(), _i32()
        self._check(lib.arm_stage_grid_coordinates(h, ctypes.byref(grid), ctypes.byref(grid_count)))
        grid_parts, grid_part_count = _i32p(), _i32()
        self._check(lib.arm_stage_grid_parts(h, ctypes.byref(grid_parts), ctypes.byref(grid_part_count)))
        uv_count = _i32()
        self._check(lib.arm_stage_uv_set_count(h, ctypes.byref(uv_count)))
        uv_sets = []
        for index in range(uv_count.value):
            name, corner_uv, uv, n = ctypes.c_char_p(), _i32p(), _f32p(), _i32()
            self._check(lib.arm_stage_uv_set(h, index, ctypes.byref(name), ctypes.byref(corner_uv),
                                             ctypes.byref(uv), ctypes.byref(n)))
            uv_sets.append((name.value.decode("utf-8"), _copy(corner_uv, view.corner_count, np.int32),
                            _copy(uv, n.value, np.float32, 2)))
        has_normals = normal_count.value > 0
        return StageData(
            stage=stage.value,
            points=_copy(view.points_xyz, view.point_count, np.float32, 3),
            face_starts=_copy(view.face_starts, view.face_count + 1, np.int32),
            corner_points=_copy(view.corner_points, view.corner_count, np.int32),
            materials=face_values(FACE_MATERIAL),
            face_groups=face_values(FACE_GROUP),
            smoothing_groups=face_values(FACE_SMOOTHING_GROUPS),
            mirror_sides=face_values(FACE_MIRROR_SIDE),
            face_parts=face_values(FACE_PART),
            part_values=_copy(parts, part_count.value, np.int32, 2),
            sharp_edges=_copy(pairs, edge_count.value, np.int32, 2),
            line_starts=starts,
            line_points=_copy(line_points, int(starts[-1]), np.int32),
            line_kinds=_copy(kinds, kind_count.value, np.int32),
            corner_normal=_copy(corner_normal, view.corner_count, np.int32) if has_normals else None,
            normals=_copy(normals, normal_count.value, np.float32, 3) if has_normals else None,
            colors=_copy(colors, color_count.value, np.uint8, 4) if color_count.value else None,
            density=_copy(density, density_count.value, np.float32) if density_count.value else None,
            uv_sets=uv_sets,
            curvature_directions=_copy(directions, direction_count.value, np.float32, 3)
            if direction_count.value else None,
            cross_directions=_copy(crosses, cross_count.value, np.float32, 3) if cross_count.value else None,
            singularities=_copy(singular, singular_count.value, np.int32, 2) if singular_count.value else None,
            cut_edges=_copy(cuts, cut_count.value, np.int32, 2) if cut_count.value else None,
            size_field=_copy(sizes, size_count.value, np.float32, 2) if size_count.value else None,
            size_limits=_copy(limits, limit_count.value, np.float32, 4) if limit_count.value else None,
            grid_uv=_copy(grid, grid_count.value, np.float32, 2) if grid_count.value else None,
            grid_parts=_copy(grid_parts, grid_part_count.value, np.int32) if grid_part_count.value else None,
        )


class CoreJob:
    """One run of the core up to a stage, in the calling thread (run) or in a thread of its own (start).

    Every call on the job's context comes from the thread that computes; another
    thread only asks it to stop (cancel), as the C interface allows. Progress,
    the end and the result are read with `fraction`, `stage`, `finished` and result().
    With a prep cache the run takes its preparation from there when the input and the
    settings the preparation reads are the same (`preparation_reused`).
    """

    def __init__(self, library: CoreLibrary, settings: CoreSettings, mesh: MeshInput, stage: int,
                 prep_cache: CorePrepCache | None = None):
        self._library = library
        self._settings = settings
        self._mesh: MeshInput | None = mesh
        self._until = stage
        self._prep_cache = prep_cache
        self.preparation_reused = False
        self._lock = threading.Lock()  # guards the context between cancel() and its creation and release
        self._core: CoreContext | None = None
        self._cancelled = False
        self._finished = threading.Event()
        self._data: StageData | None = None
        self._error: CoreError | None = None
        self.fraction = 0.0  # of the whole run, 0..1
        self.stage = 0       # arm_stage the core is working on

    @property
    def finished(self) -> bool:
        return self._finished.is_set()

    @property
    def cancelled(self) -> bool:
        return self._cancelled

    def run(self) -> None:
        """Computes in this thread; the result follows from result()."""
        self._work()

    def start(self) -> None:
        """Computes in a thread of its own."""
        threading.Thread(target=self._work, name="ArekoRemesher core", daemon=True).start()

    def wait(self, timeout: float | None = None) -> bool:
        return self._finished.wait(timeout)

    def cancel(self) -> None:
        """Asks the run to stop; safe from any thread, also before the run has begun."""
        with self._lock:
            self._cancelled = True
            if self._core is not None:
                self._core.cancel()

    def result(self) -> StageData:
        """The stage data of the finished run; raises the error that ended it."""
        if self._error is not None:
            raise self._error
        if self._data is None:
            raise CoreError("the run has not finished")
        return self._data

    def _progress(self, fraction: float, stage: int) -> None:
        self.fraction = fraction
        self.stage = stage
        if self._cancelled:
            # The core clears a request that came before arm_run began; ask again.
            self.cancel()

    def _work(self) -> None:
        try:
            with self._lock:
                if self._cancelled:
                    raise CoreRunError(ARM_CANCELLED, "cancelled")
                self._core = self._library.create_context()
            core = self._core
            if self._prep_cache is not None:
                core.use_prep_cache(self._prep_cache)
            core.set_settings(self._settings)
            core.load_input(self._mesh)
            core.run(self._until, self._progress)
            if self._cancelled:
                raise CoreRunError(ARM_CANCELLED, "cancelled")
            self.preparation_reused = core.preparation_reused
            self._data = core.stage_data()
        except CoreError as error:
            self._error = error
        except Exception as error:  # noqa: BLE001 - an error must not end the thread unseen
            traceback.print_exc()
            self._error = CoreError(str(error))
        finally:
            with self._lock:
                if self._core is not None:
                    self._core.close()
                self._core = None
            self._mesh = None
            self._prep_cache = None
            self._finished.set()


@dataclass(frozen=True)
class CoreStatus:
    """Result of checking the core once: usable or not, and why."""

    ready: bool
    version: str = ""
    reason: str = ""

    @classmethod
    def probe(cls, library: CoreLibrary | None, error: str = "") -> CoreStatus:
        if library is None:
            return cls(ready=False, reason=error)
        try:
            with library.create_context():
                pass
            return cls(ready=True, version=library.version)
        except CoreError as failure:
            return cls(ready=False, reason=str(failure))
