# Guide curves: a list of Blender curve objects per scene, each with its kind,
# on/off and width, and their points for the core (object space of the mesh).

from __future__ import annotations

import bpy
import numpy as np
from bpy.props import BoolProperty, CollectionProperty, EnumProperty, FloatProperty, IntProperty, PointerProperty
from mathutils import Vector
from mathutils.geometry import interpolate_bezier

from .kern import CURVE_BAND, CURVE_CREASE, CURVE_EDGE_LOOP, CURVE_OFF, CurveInput
from .sitzung import AddonSession
from .texte import SETTINGS_CONTEXT

# identifier, core kind, English name, English description
KINDS = (
    ("BAND", CURVE_BAND, "Band", "Quads beside the curve follow its direction"),
    ("EDGE_LOOP", CURVE_EDGE_LOOP, "Edge Loop", "An edge loop of the result follows the curve"),
    ("CREASE", CURVE_CREASE, "Crease", "An edge loop follows the curve and stays sharp"),
    ("OFF", CURVE_OFF, "Off", "The curve is ignored"),
)
_KIND_ITEMS = [(identifier, name, description, number) for identifier, number, name, description in KINDS]
_KIND_NUMBER = {identifier: number for identifier, number, _name, _description in KINDS}

# The curves decide the half mode of the symmetry: an open stage object from
# the stage "mirror half" (arm_stage 3) on is computed again.
FIRST_STAGE = 3


def _refresh(_self, context) -> None:
    session = AddonSession.current()
    if session is not None and session.stage_refresh is not None:
        session.stage_refresh.request(context, FIRST_STAGE)


def _is_curve(_self, obj) -> bool:
    return obj.type == "CURVE"


class AREKO_PG_guide_curve(bpy.types.PropertyGroup):
    """One entry of the guide curve list"""

    curve: PointerProperty(name="Curve", description="The curve object that guides the edge flow",
                           type=bpy.types.Object, poll=_is_curve, update=_refresh,
                           translation_context=SETTINGS_CONTEXT)
    kind: EnumProperty(name="Kind", description="How the curve guides the result", items=_KIND_ITEMS,
                       default="EDGE_LOOP", update=_refresh, translation_context=SETTINGS_CONTEXT)
    active: BoolProperty(name="On", description="On: the curve guides the result. Off: it is ignored", default=True,
                         update=_refresh, translation_context=SETTINGS_CONTEXT)
    automatic_width: BoolProperty(name="Automatic Width", description="The width follows from the length of the curve",
                                  default=True, update=_refresh, translation_context=SETTINGS_CONTEXT)
    width: FloatProperty(name="Width",
                         description="How far beside the curve the quads follow its direction, in object units",
                         default=0.1, min=0.0, subtype="DISTANCE", update=_refresh,
                         translation_context=SETTINGS_CONTEXT)


class AREKO_PG_guide_curves(bpy.types.PropertyGroup):
    """Guide curves of a scene"""

    use: BoolProperty(name="Use Guide Curves",
                      description="On: the curves of the list guide the result. Off: they are ignored", default=True,
                      update=_refresh, translation_context=SETTINGS_CONTEXT)
    new_kind: EnumProperty(name="Kind for New Curves", description="Kind that added curves get", items=_KIND_ITEMS,
                           default="EDGE_LOOP", translation_context=SETTINGS_CONTEXT)
    entries: CollectionProperty(type=AREKO_PG_guide_curve)
    index: IntProperty(default=0)


class AREKO_UL_guide_curves(bpy.types.UIList):
    """List of the guide curves"""

    def draw_item(self, _context, layout, _data, item, _icon, _active_data, _active_property, _index=0, _flt_flag=0):
        row = layout.row(align=True)
        row.prop(item, "active", text="")
        row.label(text=item.curve.name if item.curve is not None else "–", icon="CURVE_DATA", translate=False)


class AREKO_OT_add_guide_curves(bpy.types.Operator):
    """Adds the selected curve objects as guide curves"""

    bl_idname = "areko.add_guide_curves"
    bl_label = "Add Curves"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return any(obj.type == "CURVE" for obj in context.selected_objects)

    def execute(self, context):
        curves = context.scene.areko_curves
        present = {entry.curve for entry in curves.entries if entry.curve is not None}
        for obj in context.selected_objects:
            if obj.type != "CURVE" or obj in present:
                continue
            entry = curves.entries.add()
            entry.curve = obj
            entry.kind = curves.new_kind
            curves.index = len(curves.entries) - 1
        _refresh(None, context)
        return {"FINISHED"}


class AREKO_OT_remove_guide_curve(bpy.types.Operator):
    """Removes the chosen guide curve from the list"""

    bl_idname = "areko.remove_guide_curve"
    bl_label = "Remove Curve"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        curves = context.scene.areko_curves
        return 0 <= curves.index < len(curves.entries)

    def execute(self, context):
        curves = context.scene.areko_curves
        curves.entries.remove(curves.index)
        curves.index = max(0, min(curves.index, len(curves.entries) - 1))
        _refresh(None, context)
        return {"FINISHED"}


class GuideCurveInput:
    """Points of the guide curves of a scene in the object space of a mesh object."""

    @classmethod
    def collect(cls, context, mesh_object: bpy.types.Object) -> list[CurveInput]:
        curves = context.scene.areko_curves
        if not curves.use:
            return []
        to_mesh = mesh_object.matrix_world.inverted()
        out = []
        for entry in curves.entries:
            obj = entry.curve
            if obj is None or obj.type != "CURVE":
                continue
            matrix = to_mesh @ obj.matrix_world
            for points in cls._splines(obj.data):
                ring = np.array([tuple(matrix @ p) for p in points], dtype=np.float32)
                if len(ring) < 2:
                    continue
                out.append(CurveInput(points=ring, kind=_KIND_NUMBER[entry.kind], active=entry.active,
                                      width=-1.0 if entry.automatic_width else float(entry.width)))
        return out

    @staticmethod
    def _splines(data: bpy.types.Curve):
        """Points of every spline in curve space; a closed spline repeats its first point."""
        for spline in data.splines:
            if spline.type == "BEZIER":
                knots = spline.bezier_points
                count = len(knots)
                steps = max(spline.resolution_u, 1) + 1
                pieces = count if spline.use_cyclic_u else count - 1
                points = []
                for i in range(pieces):
                    a = knots[i]
                    b = knots[(i + 1) % count]
                    segment = interpolate_bezier(a.co, a.handle_right, b.handle_left, b.co, steps)
                    points.extend(segment if i == 0 else segment[1:])
                if count == 1:
                    points = [knots[0].co.copy()]
            else:
                points = [Vector(p.co[:3]) for p in spline.points]
                if spline.use_cyclic_u and points:
                    points.append(points[0].copy())
            if spline.use_cyclic_u and points and (points[0] - points[-1]).length > 0.0:
                points.append(points[0].copy())
            elif spline.use_cyclic_u and points:
                points[-1] = points[0].copy()
            yield points


CLASSES = (AREKO_PG_guide_curve, AREKO_PG_guide_curves, AREKO_UL_guide_curves, AREKO_OT_add_guide_curves,
           AREKO_OT_remove_guide_curve)


def register_properties() -> None:
    bpy.types.Scene.areko_curves = bpy.props.PointerProperty(type=AREKO_PG_guide_curves)


def unregister_properties() -> None:
    del bpy.types.Scene.areko_curves

