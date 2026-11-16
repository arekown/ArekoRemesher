# Stage data -> Blender object.

from __future__ import annotations

import colorsys
from dataclasses import replace

import bpy
import gpu
import numpy as np
from gpu_extras.batch import batch_for_shader

from .eingabe import FACE_SET_ATTRIBUTE
from .kern import LINE_GUIDE_CURVE, StageData
from .feldgitter import FieldGrid
from .texte import iface

# Viewport colours of the two halves of a mirrored stage (linear RGBA).
HALF_COLORS = (("Computed Half", (0.80, 0.42, 0.16, 1.0)), ("Mirrored Half", (0.18, 0.45, 0.80, 1.0)))

# Feature lines: colour of the design (#ff2d6f) and width in pixels.
LINE_COLOR = (1.0, 0.176, 0.435, 1.0)
LINE_WIDTH = 2.0
# Lines and points are drawn closer to the viewer than the faces they lie on,
# by this many pixel sizes at their own distance, so that they stay in front
# of their faces from near and far. A point needs more than a line: its disc
# faces the viewer while the face below it is tilted or lies in a dent; points
# on faces that look away from the viewer are left out instead.
LINE_DEPTH_PIXELS = 3.0
POINT_DEPTH_PIXELS = 14.0
# Guide curves cut into the mesh (#2f8cff).
GUIDE_COLOR = (0.184, 0.549, 1.0, 1.0)
# Curvature strokes on every triangle with a direction: along the direction of
# least curvature (#ffb02e) and across it (#2ec4d6), each reaching this share
# of the way from the centre of its triangle to the border, to both sides.
ALONG_COLOR = (1.0, 0.69, 0.18, 1.0)
ACROSS_COLOR = (0.18, 0.769, 0.839, 1.0)
STROKE_SHARE = 0.9
STROKE_WIDTH = 1.5
# Cross field: one cross per triangle in one colour, and the singular points by
# the turn of the crosses around them in quarter turns: 1 = three quads will
# meet there (#ff3b30), -1 = five (#2f6bff), every other turn (#bf5af2).
# Point size in pixels.
CROSS_COLOR = ALONG_COLOR
SINGULAR_COLORS = ((1, (1.0, 0.231, 0.188, 1.0)), (-1, (0.184, 0.42, 1.0, 1.0)))
OTHER_SINGULAR_COLOR = (0.749, 0.353, 0.949, 1.0)
SINGULAR_SIZE = 9.0
# Cuts: the edges along which the parts are cut open into discs (#ffe100).
CUT_COLOR = (1.0, 0.882, 0.0, 1.0)

# Part colours: hues around the colour circle in steps of the golden angle.
PART_COLORS = 24
STAGE_PARTS = 5


def part_color(i: int) -> tuple:
    hue = (i * 0.381966) % 1.0
    r, g, b = colorsys.hsv_to_rgb(hue, 0.62, 0.86)
    return (r, g, b, 1.0)


# Size field: one colour per triangle for the side length of its quads, from red (small)
# over yellow, green and cyan to blue (largest), in this many steps.
SIZE_COLORS = 32
# Custom property on the materials that colour a stage; they are removed with its object.
STAGE_COLOR_MARK = "areko_stage_color"


def size_color(i: int) -> tuple:
    hue = (2.0 / 3.0) * i / (SIZE_COLORS - 1)
    r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 0.9)
    return (r, g, b, 1.0)


# Triangles of a part that has no field to show (no size field, no grid coordinates).
NO_FIELD_COLOR = ("No Field", (0.5, 0.5, 0.5, 1.0))
# Parameterization: the grid lines (#1a3a5c) and the fields of the cells between them (#8fbddb),
# in the order of the colour numbers of the field grid; the yellow cuts stand out on both.
GRID_COLORS = ((0.102, 0.227, 0.361, 1.0), (0.561, 0.741, 0.859, 1.0))
# Raw mesh and quads: quads in a neutral grey, triangles in amber (#ffb02e), faces of five and
# more corners in the colour of the design (#ff2d6f), in the order of the face kinds below.
STAGE_RAW_MESH = 12
STAGE_QUADS = 13
FACE_KIND_COLORS = (
    ("Quad Faces", (0.62, 0.62, 0.62, 1.0)),
    ("Triangle Faces", (1.0, 0.69, 0.18, 1.0)),
    ("Faces with 5+ Corners", (1.0, 0.176, 0.435, 1.0)),
)


def _redraw_viewports() -> None:
    for window in bpy.context.window_manager.windows:
        for area in window.screen.areas:
            if area.type == "VIEW_3D":
                area.tag_redraw()


class StageLines:
    """Draws the lines of a stage over its object in every 3D viewport."""

    def __init__(self):
        self._handle = None
        self._object_name = ""
        self._shader = None
        self._point_shader = None
        self._batches = []  # (batch, colour, width)
        self._points = []  # (positions, normals, colour)

    @staticmethod
    def segments(data: StageData, guide: bool | None = None) -> np.ndarray:
        """Start and end point (object coordinates) of every line piece, as (2k, 3);
        with `guide` only the lines of guide curves (True) or the other lines (False)."""
        points = data.line_points
        if len(points) < 2:
            return np.empty((0, 3), dtype=np.float32)
        same_line = np.ones(len(points) - 1, dtype=bool)
        breaks = data.line_starts[1:-1] - 1
        same_line[breaks[(breaks >= 0) & (breaks < len(same_line))]] = False
        if guide is not None:
            count = len(data.line_starts) - 1
            kinds = data.line_kinds if len(data.line_kinds) == count else np.zeros(count, dtype=np.int32)
            line_of_piece = np.searchsorted(data.line_starts, np.arange(len(points) - 1), side="right") - 1
            same_line &= (kinds[np.clip(line_of_piece, 0, count - 1)] == LINE_GUIDE_CURVE) == guide
        pairs = np.stack((points[:-1][same_line], points[1:][same_line]), axis=1).reshape(-1)
        return data.points[pairs]

    @classmethod
    def curvature_strokes(cls, data: StageData) -> tuple[np.ndarray, np.ndarray]:
        """Strokes along and across the curvature direction of every triangle that has one, each as (2k, 3)."""
        return cls.face_strokes(data, data.curvature_directions)

    @classmethod
    def cross_strokes(cls, data: StageData) -> np.ndarray:
        """Both strokes of the cross of every triangle that has one, as (2k, 3)."""
        return np.concatenate(cls.face_strokes(data, data.cross_directions))

    @staticmethod
    def cut_segments(data: StageData) -> np.ndarray:
        """Start and end point (object coordinates) of every cut edge, as (2k, 3)."""
        if data.cut_edges is None or len(data.cut_edges) == 0:
            return np.empty((0, 3), dtype=np.float32)
        return data.points[data.cut_edges.reshape(-1)]

    @staticmethod
    def singular_points(data: StageData, turn: int | None) -> np.ndarray:
        """Numbers of the singular points with this turn among the points of the stage mesh.

        None selects the points whose turn has no colour of its own.
        """
        pairs = data.singularities
        if pairs is None or len(pairs) == 0:
            return np.empty(0, dtype=np.int32)
        if turn is None:
            wanted = ~np.isin(pairs[:, 1], [own for own, _ in SINGULAR_COLORS])
        else:
            wanted = pairs[:, 1] == turn
        return pairs[wanted, 0]

    @staticmethod
    def marks_faces(data: StageData) -> bool:
        """True when the stage has strokes on its triangles or cut edges between them."""
        return (data.curvature_directions is not None or data.cross_directions is not None
                or data.cut_edges is not None)

    @staticmethod
    def face_strokes(data: StageData, directions: np.ndarray | None) -> tuple[np.ndarray, np.ndarray]:
        """Strokes along a direction per triangle and across it, for every triangle whose direction is not zero."""
        empty = np.empty((0, 3), dtype=np.float32)
        if directions is None or len(directions) != len(data.face_starts) - 1:
            return empty, empty
        corners = data.corner_points[data.face_starts[:-1, None] + np.arange(3)]
        p = data.points[corners].astype(np.float64)
        normal = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        d = directions.astype(np.float64)
        d_len = np.linalg.norm(d, axis=1)
        n_len = np.linalg.norm(normal, axis=1)
        keep = (d_len > 0.0) & (n_len > 0.0)
        if not np.any(keep):
            return empty, empty
        p, normal, n_len = p[keep], normal[keep], n_len[keep]
        d = d[keep] / d_len[keep, None]
        across = np.cross(normal / n_len[:, None], d)
        a_len = np.linalg.norm(across, axis=1)
        across = np.divide(across, a_len[:, None], out=np.zeros_like(across), where=a_len[:, None] > 0.0)
        # A stroke from the centre leaves its triangle where the weight of one corner
        # reaches zero; `slope` is how fast the weight of each corner changes along the face.
        slope = np.cross(normal[:, None, :], np.roll(p, -2, axis=1) - np.roll(p, -1, axis=1)) \
            / (n_len * n_len)[:, None, None]
        centre = p.mean(axis=1)

        def pairs(v):
            speed = np.abs(np.einsum("fcx,fx->fc", slope, v)).max(axis=1)
            half = (STROKE_SHARE / 3.0 * np.divide(1.0, speed, out=np.zeros_like(speed), where=speed > 0.0))[:, None]
            return np.stack((centre - v * half, centre + v * half), axis=1).reshape(-1, 3).astype(np.float32)

        return pairs(d), pairs(across)

    def show(self, obj: bpy.types.Object, data: StageData) -> None:
        self.clear()
        self._shader = gpu.shader.from_builtin("POLYLINE_UNIFORM_COLOR")
        for guide, color in ((False, LINE_COLOR), (True, GUIDE_COLOR)):
            coords = self.segments(data, guide)
            if len(coords):
                self._batches.append((batch_for_shader(self._shader, "LINES", {"pos": coords}), color, LINE_WIDTH))
        for coords, color in zip(self.curvature_strokes(data), (ALONG_COLOR, ACROSS_COLOR)):
            if len(coords):
                self._batches.append((batch_for_shader(self._shader, "LINES", {"pos": coords}), color, STROKE_WIDTH))
        crosses = self.cross_strokes(data)
        if len(crosses):
            self._batches.append((batch_for_shader(self._shader, "LINES", {"pos": crosses}), CROSS_COLOR, STROKE_WIDTH))
        cuts = self.cut_segments(data)
        if len(cuts):
            self._batches.append((batch_for_shader(self._shader, "LINES", {"pos": cuts}), CUT_COLOR, LINE_WIDTH))
        for turn, color in (*SINGULAR_COLORS, (None, OTHER_SINGULAR_COLOR)):
            numbers = self.singular_points(data, turn)
            if len(numbers):
                if self._point_shader is None:
                    self._point_shader = gpu.shader.from_builtin("POINT_UNIFORM_COLOR")
                    normals = np.empty((len(obj.data.vertices), 3), dtype=np.float32)
                    obj.data.vertex_normals.foreach_get("vector", normals.reshape(-1))
                self._points.append((data.points[numbers], normals[numbers], color))
        if not self._batches and not self._points:
            self._shader = None
            return
        self._object_name = obj.name
        self._handle = bpy.types.SpaceView3D.draw_handler_add(self._draw, (), "WINDOW", "POST_VIEW")
        _redraw_viewports()

    def clear(self) -> None:
        if self._handle is not None:
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, "WINDOW")
            _redraw_viewports()
        self._handle = None
        self._object_name = ""
        self._shader = None
        self._point_shader = None
        self._batches = []
        self._points = []

    def _draw(self) -> None:
        obj = bpy.data.objects.get(self._object_name)
        region = bpy.context.region
        if obj is None or region is None or not obj.visible_get() or not (self._batches or self._points):
            return
        shader = self._shader
        depth_test, blend = gpu.state.depth_test_get(), gpu.state.blend_get()
        gpu.state.depth_test_set("LESS_EQUAL")
        gpu.state.blend_set("ALPHA")
        with gpu.matrix.push_pop(), gpu.matrix.push_pop_projection():
            gpu.matrix.multiply_matrix(obj.matrix_world)
            projection = gpu.matrix.get_projection_matrix()
            if self._batches:
                gpu.matrix.load_projection_matrix(self._closer(projection, region, LINE_DEPTH_PIXELS))
                shader.uniform_float("viewportSize", (region.width, region.height))
                for batch, color, width in self._batches:
                    shader.uniform_float("lineWidth", width)
                    shader.uniform_float("color", color)
                    batch.draw(shader)
            if self._points:
                gpu.matrix.load_projection_matrix(self._closer(projection, region, POINT_DEPTH_PIXELS))
                gpu.state.point_size_set(SINGULAR_SIZE)
                for coords, normals, color in self._points:
                    seen = coords[self._facing(coords, normals, projection)]
                    if len(seen):
                        self._point_shader.uniform_float("color", color)
                        batch_for_shader(self._point_shader, "POINTS", {"pos": seen}).draw(self._point_shader)
                gpu.state.point_size_set(1.0)
        gpu.state.blend_set(blend)
        gpu.state.depth_test_set(depth_test)

    @staticmethod
    def _facing(coords: np.ndarray, normals: np.ndarray, projection) -> np.ndarray:
        """True for every point (object coordinates) whose normal looks towards the viewer."""
        to_object = gpu.matrix.get_model_view_matrix().inverted()
        if projection[3][3] == 0.0:
            to_viewer = np.array(to_object.translation, dtype=np.float32) - coords
        else:
            to_viewer = np.array(to_object.col[2].xyz, dtype=np.float32)
        return (normals * to_viewer).sum(axis=1) > 0.0

    @staticmethod
    def _closer(projection, region, pixels: float):
        """The projection with everything `pixels` pixel sizes closer to the viewer."""
        closer = projection.copy()
        step = pixels * 2.0 / (projection[0][0] * region.width)
        if projection[3][3] == 0.0:
            # Perspective: a pixel grows with the distance, so does the step.
            closer[2][3] = projection[2][3] / (1.0 - step)
        else:
            closer[2][3] = projection[2][3] + projection[2][2] * step
        return closer


class StageObject:
    """Builds the object that shows one stage at the place of its input."""

    @staticmethod
    def _materials(entries) -> list:
        """One material per (label, colour), created where it is missing and marked as colour of a stage."""
        materials = []
        for label, color in entries:
            material = bpy.data.materials.get(label) or bpy.data.materials.new(label)
            material.diffuse_color = color
            material[STAGE_COLOR_MARK] = True
            materials.append(material)
        return materials

    @classmethod
    def _half_materials(cls) -> list:
        return cls._materials((iface(name), color) for name, color in HALF_COLORS)

    @classmethod
    def _part_materials(cls) -> list:
        return cls._materials(
            (iface("Part Color {}").format(i + 1), part_color(i)) for i in range(PART_COLORS)
        )

    @classmethod
    def _no_field_material(cls) -> tuple:
        return iface(NO_FIELD_COLOR[0]), NO_FIELD_COLOR[1]

    @classmethod
    def _size_materials(cls) -> list:
        return cls._materials(
            [(iface("Size Color {}").format(i + 1), size_color(i)) for i in range(SIZE_COLORS)]
            + [cls._no_field_material()]
        )

    @classmethod
    def _grid_materials(cls) -> list:
        """Materials in the order of the colour numbers of the field grid."""
        return cls._materials(
            [(iface("Grid Color {}").format(i + 1), color) for i, color in enumerate(GRID_COLORS)]
            + [cls._no_field_material()]
        )

    @classmethod
    def _face_kind_materials(cls) -> list:
        return cls._materials((iface(name), color) for name, color in FACE_KIND_COLORS)

    @staticmethod
    def face_kinds(data: StageData) -> np.ndarray:
        """Kind of every face: 0 quad, 1 triangle, 2 any other corner count (five and more)."""
        corners = np.diff(data.face_starts)
        return np.where(corners == 4, 0, np.where(corners == 3, 1, 2)).astype(np.int32)

    @staticmethod
    def size_steps(data: StageData) -> np.ndarray:
        """Colour step of every triangle of a size field: 0 = smallest quads, SIZE_COLORS - 1 = largest.

        The side length of a triangle's quads (the square with the area of its two sizes) is placed
        between the middle of the smallest and the largest size of its part (their geometric mean) and
        the largest, on a logarithmic scale; lengths outside take the nearer end. The lower half of the
        range is red as a whole: a field gets there in few places only, and a scale down to the smallest
        size would leave nearly all of a mesh between green and blue. A part whose two limits are equal
        is at the large end, and at the small end only where a length lies below them. The triangles of
        a part without size field (its largest size is 0) get the step SIZE_COLORS: no colour of the scale.
        """
        sizes = data.size_field.astype(np.float64)
        length = np.sqrt(np.abs(sizes[:, 0] * sizes[:, 1]))
        parts = data.face_parts if data.face_parts is not None else np.zeros(len(length), dtype=np.int32)
        low = data.size_limits[parts, 1].astype(np.float64)
        high = data.size_limits[parts, 2].astype(np.float64)
        spread = (low > 0.0) & (high > low * (1.0 + 1e-6)) & (length > 0.0)
        place = np.where(length < low * (1.0 - 1e-4), 0.0, 1.0)
        middle = np.sqrt(low[spread] * high[spread])
        place[spread] = np.clip(np.log(length[spread] / middle) / np.log(high[spread] / middle), 0.0, 1.0)
        steps = np.rint(place * (SIZE_COLORS - 1)).astype(np.int32)
        steps[~(high > 0.0)] = SIZE_COLORS
        return steps

    @staticmethod
    def _sharp_edge_mask(mesh: bpy.types.Mesh, pairs: np.ndarray) -> np.ndarray:
        edges = np.empty(len(mesh.edges) * 2, dtype=np.int64)
        mesh.edges.foreach_get("vertices", edges)
        edges = np.sort(edges.reshape(-1, 2), axis=1)
        width = np.int64(len(mesh.vertices) + 1)
        keys = edges[:, 0] * width + edges[:, 1]
        wanted = np.sort(pairs.astype(np.int64), axis=1)
        wanted_keys = np.sort(wanted[:, 0] * width + wanted[:, 1])
        if len(wanted_keys) == 0:
            return np.zeros(len(keys), dtype=bool)
        pos = np.clip(np.searchsorted(wanted_keys, keys), 0, len(wanted_keys) - 1)
        return wanted_keys[pos] == keys

    @classmethod
    def build_mesh(cls, name: str, data: StageData, materials: list, smooth: bool = True) -> bpy.types.Mesh:
        mesh = bpy.data.meshes.new(name)
        mesh.vertices.add(len(data.points))
        mesh.vertices.foreach_set("co", data.points.reshape(-1))
        mesh.loops.add(len(data.corner_points))
        mesh.loops.foreach_set("vertex_index", data.corner_points)
        mesh.polygons.add(len(data.face_starts) - 1)
        mesh.polygons.foreach_set("loop_start", data.face_starts[:-1])
        mesh.update(calc_edges=True)
        for material in materials:
            mesh.materials.append(material)
        if data.materials is not None:
            mesh.polygons.foreach_set("material_index", data.materials)
        if data.colors is not None:
            colors = mesh.color_attributes.new("Color", "BYTE_COLOR", "POINT")
            colors.data.foreach_set("color_srgb", (data.colors.astype(np.float32) / 255.0).reshape(-1))
        if smooth:
            mesh.shade_smooth()
        else:
            mesh.shade_flat()
        if data.normals is not None:
            # Faces whose corners share one normal are flat; the rest takes the normals per corner.
            index = data.corner_normal
            starts = data.face_starts[:-1]
            first = np.repeat(index[starts], np.diff(data.face_starts))
            mixed = np.add.reduceat((index != first).astype(np.int32), starts) > 0
            flat = mesh.attributes.new("sharp_face", "BOOLEAN", "FACE")
            flat.data.foreach_set("value", ~mixed)
            mesh.normals_split_custom_set(data.normals[index])
        # Setting custom normals marks edges sharp on its own; the stage keeps exactly its sharp edges.
        mask = cls._sharp_edge_mask(mesh, data.sharp_edges)
        if mask.any() or "sharp_edge" in mesh.attributes:
            sharp = mesh.attributes.get("sharp_edge") or mesh.attributes.new("sharp_edge", "BOOLEAN", "EDGE")
            sharp.data.foreach_set("value", mask)
        mesh.update()
        return mesh

    @classmethod
    def create(cls, context, source: bpy.types.Object, data: StageData, name: str) -> bpy.types.Object:
        smooth = True
        if data.mirror_sides is not None:
            # A mirrored stage shows its two halves in two colours instead of the input's materials.
            data.materials = data.mirror_sides
            materials = cls._half_materials()
        elif data.stage == STAGE_PARTS and data.face_parts is not None:
            # The stage "parts" shows every part in its own colour.
            data.materials = (data.face_parts % PART_COLORS).astype(np.int32)
            materials = cls._part_materials()
        elif data.size_field is not None and data.size_limits is not None:
            # The stage "size field" shows the size of the quads as a colour per triangle.
            data.materials = cls.size_steps(data)
            materials = cls._size_materials()
        elif data.grid_uv is not None:
            # The stage "parameterization" shows the cells of the grid as fields between its
            # lines: its triangles are cut along the grid lines, every piece takes the colour of
            # a field or of a line.
            data = FieldGrid.build(data)
            materials = cls._grid_materials()
        elif data.stage in (STAGE_RAW_MESH, STAGE_QUADS):
            # The stages "raw mesh" and "quads" show their faces by corner count: quads stay
            # neutral, triangles and faces of five and more corners stand out. Every face is shaded
            # flat in its own plane: a degenerate face would spoil the smooth normals of its corners
            # and turn its sound neighbours black, and smoothing would hide folds of the mesh.
            data.materials = cls.face_kinds(data)
            materials = cls._face_kind_materials()
            smooth = False
        else:
            materials = [slot.material for slot in source.material_slots]
        obj = bpy.data.objects.new(name, cls.build_mesh(name, data, materials, smooth))
        cls.place_like(context, source, obj)
        return obj

    @staticmethod
    def place_like(context, source: bpy.types.Object, obj: bpy.types.Object) -> None:
        """Links `obj` into the first collection of `source` and puts it at the same place in the world."""
        collections = source.users_collection or (context.scene.collection,)
        collections[0].objects.link(obj)
        obj.matrix_world = source.matrix_world.copy()

    @staticmethod
    def remove(obj: bpy.types.Object) -> None:
        mesh = obj.data
        colors = [] if mesh is None else [m for m in mesh.materials if m is not None and m.get(STAGE_COLOR_MARK)]
        bpy.data.objects.remove(obj)
        if mesh is not None and mesh.users == 0:
            bpy.data.meshes.remove(mesh)
        # The colours of a stage leave the file with the last mesh that shows them.
        for material in colors:
            if material.users == 0:
                bpy.data.materials.remove(material)


# Face attributes of the result: the polygroups as Blender's face sets, the smoothing groups as bit masks.
SMOOTHING_GROUPS_ATTRIBUTE = "smoothing_groups"
# Object settings taken over from the input: locks, rotation mode and the mirror options of mesh editing.
TAKEN_OBJECT_SETTINGS = ("rotation_mode", "lock_location", "lock_rotation", "lock_rotation_w", "lock_rotations_4d",
                         "lock_scale", "use_mesh_mirror_x", "use_mesh_mirror_y", "use_mesh_mirror_z")


class ResultObject:
    """Builds the finished quad mesh (stage "result") as a new object at the place of its input."""

    @staticmethod
    def input_is_flat(source: bpy.types.Object) -> bool:
        """The input counts as flat shaded when its first face is (one value for the whole mesh)."""
        faces = source.data.polygons
        return len(faces) > 0 and not faces[0].use_smooth

    @classmethod
    def create(cls, context, source: bpy.types.Object, data: StageData, name: str,
               edge_crease: bool) -> bpy.types.Object:
        """The result object: the input's material slots with the core's material per face, face sets,
        smoothing groups, sharp edges (with `edge_crease` also as full edge crease), point colours; flat
        when the input is flat, else smooth. The result carries no normals of its own: its shading comes
        from the sharp edges."""
        data = replace(data, corner_normal=None, normals=None)
        materials = [slot.material for slot in source.material_slots]
        mesh = StageObject.build_mesh(name, data, materials, smooth=not cls.input_is_flat(source))
        cls._face_values(mesh, data)
        if edge_crease:
            cls._edge_crease(mesh)
        obj = bpy.data.objects.new(name, mesh)
        # The rotation mode comes first: placing the object splits its matrix into rotations of that mode.
        for setting in TAKEN_OBJECT_SETTINGS:
            setattr(obj, setting, getattr(source, setting))
        StageObject.place_like(context, source, obj)
        return obj

    @staticmethod
    def _face_values(mesh: bpy.types.Mesh, data: StageData) -> None:
        for attribute, values in ((FACE_SET_ATTRIBUTE, data.face_groups),
                                  (SMOOTHING_GROUPS_ATTRIBUTE, data.smoothing_groups)):
            if values is not None:
                mesh.attributes.new(attribute, "INT", "FACE").data.foreach_set("value", values)

    @staticmethod
    def _edge_crease(mesh: bpy.types.Mesh) -> None:
        """Full edge crease on every sharp edge of the result."""
        sharp = mesh.attributes.get("sharp_edge")
        if sharp is None:
            return
        flags = np.empty(len(mesh.edges), dtype=bool)
        sharp.data.foreach_get("value", flags)
        if flags.any():
            crease = mesh.attributes.new("crease_edge", "FLOAT", "EDGE")
            crease.data.foreach_set("value", flags.astype(np.float32))
