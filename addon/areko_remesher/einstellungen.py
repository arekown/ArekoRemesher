# Settings of a run, stored per scene. Field names and types follow the core's
# settings block; the defaults come from the core itself (one source). A few
# settings belong to the add-on alone (ADDON_FIELDS) and never reach the core.

from __future__ import annotations

import math
import struct

import bpy
from bpy.props import BoolProperty, EnumProperty, FloatProperty, IntProperty

from .kern import CoreSettings
from .sitzung import AddonSession
from .texte import SETTINGS_CONTEXT

# field, kind, English name, limits in core units, first stage that reads it
# (arm_stage number; None: no stage of this version reads it yet)
FIELDS = (
    ("quad_count", "INT", "Quads (Count)", {"min": 1}, 1),
    ("quad_percentage", "FLOAT", "Quads (Share of Input, %)", {"min": 0.0, "max": 1000.0}, 1),
    ("edge_length", "FLOAT", "Edge Length", {"min": 0.0}, 1),
    ("exact_count", "BOOL", "Keep Count Strictly", {}, 10),
    ("color_density", "BOOL", "Density from Color", {}, 0),
    ("attribute_density", "BOOL", "Density from Attribute", {}, 10),
    ("curvature", "FLOAT", "Curvature Densifies", {"min": 0.0, "max": 100.0}, 10),
    ("points_per_circle", "FLOAT", "Points per Circle", {"min": 0.0}, 10),
    ("density_ratio", "FLOAT", "Density Range", {"min": 0.0}, 10),
    ("max_quad_ratio", "FLOAT", "Largest Side Ratio", {"min": 1.0}, 10),
    ("follow_materials", "BOOL", "Follow Material Borders", {}, 1),
    ("follow_face_groups", "BOOL", "Follow Face Groups", {}, 1),
    ("follow_smoothing_groups", "BOOL", "Follow Smoothing Groups", {}, 1),
    ("follow_sharp_edges", "BOOL", "Follow Sharp Edges", {}, 1),
    ("follow_normal_splits", "BOOL", "Follow Normal Splits", {}, 1),
    ("detect_creases", "BOOL", "Detect Creases", {}, 4),
    ("merge_cad_groups", "BOOL", "Merge CAD Faces", {}, 1),
    ("weld_seams", "BOOL", "Weld Open Seams", {}, 1),
    ("weld_as_crease", "BOOL", "Welded Seams as Crease", {}, 1),
    ("denoise", "BOOL", "Denoise", {}, 2),
    ("crease_angle", "ANGLE", "Crease Angle", {"min": 0.0, "max": 180.0}, 2),
    ("follow_borders", "BOOL", "Follow Borders", {}, 7),
    ("freeze_borders", "BOOL", "Freeze Borders", {}, 1),
    ("border_density_factor", "FLOAT", "Border Density Factor", {"min": 0.0}, 7),
    ("mirror_x", "BOOL", "Mirror X", {}, 3),
    ("mirror_y", "BOOL", "Mirror Y", {}, 3),
    ("mirror_z", "BOOL", "Mirror Z", {}, 3),
    ("mirror_at_center", "BOOL", "Mirror at Object Center", {}, 3),
    ("detect_symmetry", "BOOL", "Find Symmetry", {}, 2),
    ("freeze_half", "BOOL", "Keep Half", {}, 3),
    ("split_at_creases", "BOOL", "Split at Creases", {}, 14),
    ("creases_to_edge_crease", "BOOL", "Creases as Edge Crease", {}, None),
    ("clean_output", "BOOL", "Clean Output", {}, 14),
    ("keep_colors", "BOOL", "Keep Colors", {}, 0),
    ("thread_count", "INT", "Cores", {"min": 0}, None),
)

# Tooltip of every field (English; German in texte.py).
DESCRIPTIONS = {
    "quad_count": "About how many quads the result gets",
    "quad_percentage": "Number of quads as a share of the input's faces",
    "edge_length": "Length of the quad edges, in object units",
    "exact_count": "Off: the number of quads may change to suit the shape. On: it is kept",
    "color_density": "The first colour attribute sets the size: red gives smaller quads, cyan larger ones",
    "attribute_density": "The point attribute \"density\" sets the size: 4 gives smaller quads, 0.25 larger ones",
    "curvature": "How much smaller the quads get where the surface bends; 0 = the same size everywhere",
    "points_per_circle": "Quads around a full bend; 0 = from Curvature Densifies",
    "density_ratio": "Ratio of the largest to the smallest quads; 0 = from Curvature Densifies, 1 = one size",
    "max_quad_ratio": "Largest ratio of the long to the short side of a quad",
    "follow_materials": "Edge loops follow the borders between materials",
    "follow_face_groups": "Edge loops follow the borders between face sets",
    "follow_smoothing_groups": "Edge loops follow the borders between smoothing groups",
    "follow_sharp_edges": "Edge loops follow the edges marked sharp",
    "follow_normal_splits": "Edge loops follow the edges where the normals split",
    "detect_creases": "Finds creases where the faces meet at a sharp angle",
    "merge_cad_groups": "Merges narrow strips between face groups, as CAD meshes have them",
    "weld_seams": "Welds open seams whose border points lie exactly on top of each other",
    "weld_as_crease": "Welded seams count as creases",
    "denoise": "Smooths out noise of the surface before the computation",
    "crease_angle": "Faces meeting at a sharper angle than this form a crease",
    "follow_borders": "Edge loops follow the open borders",
    "freeze_borders": "The points of the open borders stay where they are",
    "border_density_factor": "Scales the size of the quads along frozen borders; 1 = as elsewhere",
    "mirror_x": "The result is mirrored across the X plane",
    "mirror_y": "The result is mirrored across the Y plane",
    "mirror_z": "The result is mirrored across the Z plane",
    "mirror_at_center": "The mirror plane runs through the middle of the object; off: through its origin",
    "detect_symmetry": "Finds the symmetry of every part and computes one half",
    "freeze_half": "Holds the border at the mirror plane in place; only useful with Freeze Borders",
    "split_at_creases": "Splits the points of the result along its creases",
    "creases_to_edge_crease": "The creases of the result get full edge crease",
    "clean_output": "Dissolves leftover edges and separates points where the result is not manifold",
    "keep_colors": "The colour attribute of the input is carried over to the result",
    "thread_count": "Cores that compute the separate pieces of the mesh side by side; 0 = all. "
                    "The result is the same for every number",
    "hide_input": "Hides the input object once the result is there",
    "remember_preparation": "Keeps the prepared mesh for the next run: if only settings that act after the "
                            "preparation change, it is skipped. The result is the same",
    "size_kind": "What sets the size of the quads",
    "density_paint": "Colour of the density brush: 4 (red) gives smaller quads, 0.25 (cyan) larger ones",
}

# Exactly one of these sets the size of the quads: kind, English name, field, tooltip.
SIZE_KINDS = (
    ("COUNT", "Count", "quad_count", "The result gets about this many quads"),
    ("SHARE", "Share (%)", "quad_percentage", "Number of quads as a share of the input's faces"),
    ("LENGTH", "Edge Length", "edge_length", "Length of the quad edges, in object units"),
)
SIZE_FIELD = {kind: name for kind, _label, name, _tip in SIZE_KINDS}
# Range of the density brush (factor of the quad density, as the core reads the colours).
DENSITY_FINEST = 4.0
DENSITY_COARSEST = 0.25


def _paint_color(self, context) -> None:
    from .bedienung import DensityBrush  # the brush lives with the operators
    DensityBrush.apply(context, self.density_paint)


# Settings of the add-on alone (not in the core block): field, property, options.
ADDON_FIELDS = (
    ("hide_input", BoolProperty, {"name": "Hide Input", "default": True}),
    ("remember_preparation", BoolProperty, {"name": "Remember Preparation", "default": True}),
    ("size_kind", EnumProperty, {"name": "Size by", "default": "COUNT",
                                 "items": [(kind, label, tip) for kind, label, _name, tip in SIZE_KINDS]}),
    ("density_paint", FloatProperty, {"name": "Brush Density", "default": 1.0, "min": DENSITY_COARSEST,
                                      "max": DENSITY_FINEST, "update": _paint_color}),
)

_PROPERTY = {"INT": IntProperty, "FLOAT": FloatProperty, "ANGLE": FloatProperty, "BOOL": BoolProperty}
# extra property options per kind: an angle is shown in degrees and stored in radians
_OPTIONS = {"ANGLE": {"subtype": "ANGLE"}}


def _degrees(radians: float) -> float:
    """Degrees for the core: the shortest decimal whose radians are stored as `radians`."""
    exact = math.degrees(radians)
    stored = _single(radians)
    for digits in range(7):
        candidate = round(exact, digits)
        if _single(math.radians(candidate)) == stored:
            return candidate
    return exact


def _single(value: float) -> float:
    """`value` rounded to the single precision of a Blender float property."""
    return struct.unpack("f", struct.pack("f", value))[0]


# panel value -> core value and back; the core takes angles in degrees
_TO_CORE = {"INT": int, "FLOAT": float, "ANGLE": _degrees, "BOOL": bool}
_FROM_CORE = {"INT": int, "FLOAT": float, "ANGLE": math.radians, "BOOL": bool}


def _refresh_from(first_stage: int):
    """Update callback: an open stage object from `first_stage` on is computed again, a moment later."""
    def update(_self, context) -> None:
        session = AddonSession.current()
        if session is not None and session.stage_refresh is not None:
            session.stage_refresh.request(context, first_stage)
    return update


def _to_core(self, defaults: CoreSettings) -> CoreSettings:
    """The core block for these settings (struct_size from `defaults`)."""
    block = CoreSettings.from_buffer_copy(defaults)
    for name, kind, _label, _limits, _stage in FIELDS:
        setattr(block, name, _TO_CORE[kind](getattr(self, name)))
    # Only the chosen size kind reaches the core; the count stays as the core's fallback.
    if self.size_kind != "SHARE":
        block.quad_percentage = 0.0
    if self.size_kind != "LENGTH":
        block.edge_length = 0.0
    return block


def addon_field_names() -> tuple[str, ...]:
    return tuple(name for name, _prop, _options in ADDON_FIELDS)


def settings_group(defaults: CoreSettings | None) -> type:
    """PropertyGroup class with one property per field of the core block."""
    annotations = {}
    for name, kind, label, limits, first_stage in FIELDS:
        options = {key: _FROM_CORE[kind](value) for key, value in limits.items()}
        options.update(_OPTIONS.get(kind, {}))
        options["description"] = DESCRIPTIONS[name]
        if first_stage is not None:
            options["update"] = _refresh_from(first_stage)
        if defaults is not None:
            options["default"] = _FROM_CORE[kind](getattr(defaults, name))
            if "min" in options and options["default"] < options["min"]:
                options["min"] = options["default"]
        annotations[name] = _PROPERTY[kind](name=label, translation_context=SETTINGS_CONTEXT, **options)
    for name, prop, options in ADDON_FIELDS:
        options = dict(options, description=DESCRIPTIONS[name])
        if name == "size_kind":
            options["update"] = _refresh_from(1)
        annotations[name] = prop(translation_context=SETTINGS_CONTEXT, **options)
    return type("AREKO_PG_settings", (bpy.types.PropertyGroup,),
                {"__annotations__": annotations, "to_core": _to_core})
