# Operators of the add-on.

from __future__ import annotations

import time
from dataclasses import dataclass

import bpy
import numpy as np
from mathutils import Color

from .ausgabe import ResultObject, part_color
from .einblick import STAGES, InsightView
from .eingabe import FACE_SET_ATTRIBUTE
from .einstellungen import DENSITY_COARSEST, DENSITY_FINEST, FIELDS, addon_field_names
from .kern import ARM_CANCELLED, CoreError, CoreRunError
from .sitzung import AddonSession
from .texte import STAGE_CONTEXT, grouped, iface, report_text

# arm_stage of the finished quad mesh, arm_status of a missing mesh and of a result without faces.
STAGE_RESULT = 14
NO_MESH = 4
NO_RESULT = 5
# Seconds between two looks at a core that computes in the background.
POLL_INTERVAL = 0.1
_STAGE_NAMES = {number: label for _identifier, number, label in STAGES}


def redraw_sidebars() -> None:
    for window in bpy.context.window_manager.windows:
        for area in window.screen.areas:
            if area.type == "VIEW_3D":
                for region in area.regions:
                    if region.type == "UI":
                        region.tag_redraw()


def selected_mesh(operator: bpy.types.Operator, context) -> bpy.types.Object | None:
    """The one selected mesh object; None after an error report."""
    selected = context.selected_objects
    if len(selected) != 1:
        operator.report({"ERROR"}, report_text("Select exactly one object"))
        return None
    if selected[0].type != "MESH":
        operator.report({"ERROR"}, report_text("The selected object is not a mesh"))
        return None
    return selected[0]


def flat_with_normal_splits(context, obj: bpy.types.Object | None) -> bool:
    """Following normal splits on a flat shaded input usually gives a poor result."""
    return (obj is not None and obj.type == "MESH" and context.scene.areko_remesher.follow_normal_splits
            and ResultObject.input_is_flat(obj))


@dataclass(frozen=True)
class RunSummary:
    """What the panel shows of the last finished run."""

    name: str
    quads: int
    triangles: int
    seconds: float


class QuadRun:
    """One run of "Create Quads": the input is read at once, the core computes in this thread
    or in the background, and finish() builds the result object (nothing is created before)."""

    def __init__(self, context, source: bpy.types.Object, session: AddonSession):
        settings = context.scene.areko_remesher
        self.source_name = source.name
        self.name = iface("{} Quads").format(source.name)
        self.edge_crease = settings.creases_to_edge_crease
        self.hide_input = settings.hide_input
        self.started = time.perf_counter()
        self.job = InsightView.core_job(context, source, session, STAGE_RESULT)

    @classmethod
    def begin(cls, operator: bpy.types.Operator, context) -> QuadRun | None:
        """Checks the selection and reads the input in object mode; None after an error report."""
        # An open stage object gives its place back to the input first.
        insight = InsightView(context)
        if insight.is_open:
            insight.close()
        source = selected_mesh(operator, context)
        if source is None:
            return None
        if flat_with_normal_splits(context, source):
            operator.report({"WARNING"}, report_text("Flat shading with Follow Normal Splits usually gives a poor "
                                                     "result"))
        # The input is read in object mode (edit changes reach the mesh); the mode returns at once.
        mode = source.mode
        if mode != "OBJECT":
            bpy.ops.object.mode_set(mode="OBJECT")
        try:
            return cls(context, source, AddonSession.current())
        except CoreError as error:
            report_failure(operator, error)
            return None
        finally:
            if mode != "OBJECT":
                bpy.ops.object.mode_set(mode=mode)

    def progress_text(self) -> str:
        stage = iface(_STAGE_NAMES.get(self.job.stage, ""), STAGE_CONTEXT)
        return f"{stage} · {int(self.job.fraction * 100.0)} %"

    def finish(self, context) -> RunSummary:
        """Builds the result object of the finished job; raises the error that ended the job."""
        data = self.job.result()
        if len(data.face_starts) < 2:
            raise CoreRunError(NO_RESULT, "the result has no faces")
        source = bpy.data.objects.get(self.source_name)
        if source is None or source.type != "MESH":
            raise CoreRunError(NO_MESH, "the input object is gone")
        active = context.view_layer.objects.active
        mode = active.mode if active is not None else "OBJECT"
        if mode != "OBJECT":
            bpy.ops.object.mode_set(mode="OBJECT")
        result = ResultObject.create(context, source, data, self.name, edge_crease=self.edge_crease)
        if self.hide_input:
            source.hide_set(True)
            source.hide_viewport = True
            source.hide_render = True
        for obj in context.selected_objects:
            obj.select_set(False)
        result.select_set(True)
        context.view_layer.objects.active = result
        if mode != "OBJECT":
            bpy.ops.object.mode_set(mode=mode)
        corners = np.diff(data.face_starts)
        return RunSummary(result.name, int((corners == 4).sum()), int((corners == 3).sum()),
                          time.perf_counter() - self.started)


def report_failure(operator: bpy.types.Operator, error: CoreError) -> None:
    if isinstance(error, CoreRunError) and error.status == ARM_CANCELLED:
        operator.report({"INFO"}, report_text(error.message))
    elif isinstance(error, CoreRunError):
        operator.report({"ERROR"}, report_text(error.message))
    else:
        operator.report({"ERROR"}, f"{report_text('Internal error')}: {error}")


class AREKO_OT_create_quads(bpy.types.Operator):
    """Create a quad mesh from the selected object"""

    bl_idname = "areko.create_quads"
    bl_label = "Create Quads"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        session = AddonSession.current()
        return session is not None and session.core.ready and session.run is None

    def execute(self, context):
        # From a script or the redo panel: computes at once.
        run = QuadRun.begin(self, context)
        if run is None:
            return {"CANCELLED"}
        run.job.run()
        return self._finish(context, run)

    def invoke(self, context, event):
        # From the button or the shortcut: computes in the background, Blender stays usable.
        session = AddonSession.current()
        run = QuadRun.begin(self, context)
        if run is None:
            return {"CANCELLED"}
        session.run = run
        session.last_result = None
        run.job.start()
        self._timer = context.window_manager.event_timer_add(POLL_INTERVAL, window=context.window)
        context.window_manager.modal_handler_add(self)
        redraw_sidebars()
        return {"RUNNING_MODAL"}

    def modal(self, context, event):
        session = AddonSession.current()
        run = session.run if session is not None else None
        if run is None:
            self._end(context)
            return {"CANCELLED"}
        if event.type == "ESC" and event.value == "PRESS":
            run.job.cancel()
            return {"RUNNING_MODAL"}
        if event.type == "TIMER":
            redraw_sidebars()
            if run.job.finished:
                self._end(context)
                return self._finish(context, run)
        return {"PASS_THROUGH"}

    def cancel(self, context):
        # Blender ends the operator (e.g. a file is loaded): the run stops and leaves nothing.
        session = AddonSession.current()
        if session is not None and session.run is not None:
            session.run.job.cancel()
        self._end(context)

    def _end(self, context) -> None:
        timer = getattr(self, "_timer", None)
        if timer is not None:
            context.window_manager.event_timer_remove(timer)
            self._timer = None
        session = AddonSession.current()
        if session is not None:
            session.run = None
        redraw_sidebars()

    def _finish(self, context, run: QuadRun) -> set:
        try:
            summary = run.finish(context)
        except CoreError as error:
            report_failure(self, error)
            return {"CANCELLED"}
        session = AddonSession.current()
        if session is not None:
            session.last_result = summary
        self.report({"INFO"}, report_text("Done: {} quads").format(grouped(summary.quads)))
        return {"FINISHED"}


class AREKO_OT_cancel_quads(bpy.types.Operator):
    """Stops the running computation; nothing is created"""

    bl_idname = "areko.cancel_quads"
    bl_label = "Cancel (Esc)"

    @classmethod
    def poll(cls, context):
        session = AddonSession.current()
        return session is not None and session.run is not None

    def execute(self, context):
        AddonSession.current().run.job.cancel()
        return {"FINISHED"}


class AREKO_OT_reset_settings(bpy.types.Operator):
    """Sets all settings back to their defaults"""

    bl_idname = "areko.reset_settings"
    bl_label = "Reset"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        settings = context.scene.areko_remesher
        for name in (*(field[0] for field in FIELDS), *addon_field_names()):
            settings.property_unset(name)
        curves = context.scene.areko_curves
        for name in ("use", "new_kind"):
            curves.property_unset(name)
        self.report({"INFO"}, report_text("Settings reset"))
        return {"FINISHED"}


class AREKO_OT_face_groups_to_materials(bpy.types.Operator):
    """Gives every face set of the selected object a new material in a colour of its own"""

    bl_idname = "areko.face_groups_to_materials"
    bl_label = "Face Groups to Materials"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        obj = selected_mesh(self, context)
        if obj is None:
            return {"CANCELLED"}
        mode = obj.mode
        if mode != "OBJECT":
            bpy.ops.object.mode_set(mode="OBJECT")
        try:
            count = self._assign(obj.data)
        finally:
            if mode != "OBJECT":
                bpy.ops.object.mode_set(mode=mode)
        if count == 0:
            self.report({"ERROR"}, report_text("No face groups found"))
            return {"CANCELLED"}
        self.report({"INFO"}, report_text("{} materials created").format(count))
        return {"FINISHED"}

    @staticmethod
    def _assign(mesh: bpy.types.Mesh) -> int:
        """One new material per face set behind the existing slots; 0 when the mesh has fewer than two."""
        groups = mesh.attributes.get(FACE_SET_ATTRIBUTE)
        if groups is None or groups.domain != "FACE" or groups.data_type != "INT":
            return 0
        values = np.empty(len(mesh.polygons), dtype=np.int32)
        groups.data.foreach_get("value", values)
        distinct, index = np.unique(values, return_inverse=True)
        if len(distinct) < 2:
            return 0
        first = len(mesh.materials)
        for i, group in enumerate(distinct):
            material = bpy.data.materials.new(iface("Face Group {}").format(int(group)))
            color = part_color(i)
            material.diffuse_color = color
            if material.node_tree is not None:
                shader = next((n for n in material.node_tree.nodes if n.type == "BSDF_PRINCIPLED"), None)
                if shader is not None:
                    shader.inputs["Base Color"].default_value = color
            mesh.materials.append(material)
        mesh.polygons.foreach_set("material_index", (index + first).astype(np.int32))
        mesh.update()
        return len(distinct)


class DensityBrush:
    """The paint colour of a density factor, as the core reads the colours back."""

    @staticmethod
    def color(value: float) -> tuple[float, float, float]:
        """sRGB colour of a factor 0.25..4: white at 1, towards red (smaller quads) above, towards cyan below."""
        share = 0.0
        if value > 1.0:
            share = (value - 1.0) / (DENSITY_FINEST - 1.0)
        elif value < 1.0:
            share = -((1.0 / value) - 1.0) / ((1.0 / DENSITY_COARSEST) - 1.0)
        share = min(1.0, max(-1.0, share))
        if share > 0.0:
            return (1.0, 1.0 - share, 1.0 - share)
        if share < 0.0:
            return (1.0 + share, 1.0, 1.0)
        return (1.0, 1.0, 1.0)

    @classmethod
    def apply(cls, context, value: float) -> None:
        """Sets the colour for vertex paint and sculpt painting (brush colours are scene linear)."""
        linear = Color(cls.color(value)).from_srgb_to_scene_linear()
        tools = context.tool_settings
        for paint in (tools.vertex_paint, tools.sculpt):
            if paint is None:
                continue
            paint.unified_paint_settings.color = linear
            if paint.brush is not None and not paint.unified_paint_settings.use_unified_color:
                paint.brush.color = linear


class AREKO_OT_paint_density(bpy.types.Operator):
    """Vertex paint with the colour of Brush Density: red gives smaller quads, cyan larger ones"""

    bl_idname = "areko.paint_density"
    bl_label = "Paint Density"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == "MESH"

    def execute(self, context):
        obj = context.active_object
        if obj.mode != "OBJECT":
            bpy.ops.object.mode_set(mode="OBJECT")
        colors = obj.data.color_attributes
        if not len(colors):
            # A white map: factor 1 everywhere until painted.
            layer = colors.new(iface("Density"), "BYTE_COLOR", "POINT")
            layer.data.foreach_set("color", np.ones(len(layer.data) * 4, dtype=np.float32))
        # The core reads the first colour attribute; painting goes there.
        colors.active_color = colors[0]
        settings = context.scene.areko_remesher
        settings.color_density = True
        bpy.ops.object.mode_set(mode="VERTEX_PAINT")
        DensityBrush.apply(context, settings.density_paint)
        return {"FINISHED"}


CLASSES = (AREKO_OT_create_quads, AREKO_OT_cancel_quads, AREKO_OT_reset_settings, AREKO_OT_face_groups_to_materials,
           AREKO_OT_paint_density)
