# Insight: shows one stage of the computation as its own object in place of
# the input; the input is hidden meanwhile (never two meshes on top of each other).

from __future__ import annotations

import bpy
from bpy.props import BoolProperty, CollectionProperty, EnumProperty, IntProperty, PointerProperty, StringProperty

from .ausgabe import StageLines, StageObject
from .eingabe import BlenderMeshInput
from .kern import CoreError, CoreJob, CoreRunError
from .kurven import GuideCurveInput
from .sitzung import AddonSession
from .texte import STAGE_CONTEXT, iface, report_text

# All stages of the core (arm_stage): identifier, number, English name.
STAGES = (
    ("INPUT", 0, "Input"),
    ("PREPARED", 1, "Prepared"),
    ("CLEANED", 2, "Cleaned"),
    ("MIRROR_HALF", 3, "Mirror Half"),
    ("FEATURE_LINES", 4, "Feature Lines"),
    ("PARTS", 5, "Parts"),
    ("CURVES", 6, "Curves"),
    ("CURVATURE", 7, "Curvature"),
    ("CROSS_FIELD", 8, "Cross Field"),
    ("CUTS", 9, "Cuts"),
    ("SIZE_FIELD", 10, "Size Field"),
    ("PARAMETERIZATION", 11, "Parameterization"),
    ("RAW_MESH", 12, "Raw Mesh"),
    ("QUADS", 13, "Quads"),
    ("RESULT", 14, "Result"),
)
_BY_ID = {identifier: (number, label) for identifier, number, label in STAGES}
NO_STAGE = "–"

_items_cache: list = []


def stage_items(_self, _context):
    """Stages the loaded core can compute."""
    session = AddonSession.current()
    last = session.library.last_stage if session and session.library else -1
    _items_cache[:] = [(identifier, label, "") for identifier, number, label in STAGES if number <= last]
    return _items_cache


def stage_label(identifier: str) -> str:
    return _BY_ID[identifier][1] if identifier in _BY_ID else NO_STAGE


class AREKO_PG_stage_part(bpy.types.PropertyGroup):
    """One part of the shown stage."""

    budget: IntProperty(default=0)
    skipped: BoolProperty(default=False)


class AREKO_PG_insight(bpy.types.PropertyGroup):
    """What the insight currently shows (per scene)."""

    stage: StringProperty(default="")
    source: PointerProperty(type=bpy.types.Object)
    view: PointerProperty(type=bpy.types.Object)
    point_count: IntProperty(default=0)
    parts: CollectionProperty(type=AREKO_PG_stage_part)


class InsightView:
    """Shows, replaces and closes the stage object of a scene."""

    def __init__(self, context):
        self.context = context
        self.state = context.scene.areko_insight

    @property
    def is_open(self) -> bool:
        return self.state.view is not None

    def source(self) -> bpy.types.Object | None:
        if self.is_open:
            return self.state.source
        obj = self.context.active_object
        return obj if obj is not None and obj.type == "MESH" else None

    @staticmethod
    def core_job(context, source: bpy.types.Object, session: AddonSession, stage: int) -> CoreJob:
        """The run of the core on the evaluated mesh and guide curves of `source` up to `stage`.

        The input and the settings are read here, from Blender; the job computes later,
        in this thread (run) or in a thread of its own (start).
        """
        mesh_input = BlenderMeshInput.from_object(source, context.evaluated_depsgraph_get())
        mesh_input.curves = GuideCurveInput.collect(context, source)
        scene_settings = context.scene.areko_remesher
        settings = scene_settings.to_core(session.library.default_settings())
        return CoreJob(session.library, settings, mesh_input, stage,
                       session.prep_cache(scene_settings.remember_preparation))

    @classmethod
    def compute(cls, context, source: bpy.types.Object, session: AddonSession, stage: int):
        """Runs the core on the evaluated mesh and guide curves of `source` up to `stage`, in this thread."""
        job = cls.core_job(context, source, session, stage)
        job.run()
        return job.result()

    def show(self, identifier: str, session: AddonSession) -> bpy.types.Object:
        context = self.context
        source = self.source()
        if source is None:
            raise CoreRunError(4, "no mesh object")
        number, label = _BY_ID[identifier]
        if self.is_open:
            session.stage_lines.clear()
            StageObject.remove(self.state.view)
            self.state.view = None
            source.hide_set(False)
        data = self.compute(context, source, session, number)
        view = StageObject.create(context, source, data, iface("Stage: {}").format(iface(label, STAGE_CONTEXT)))
        # Computed stages show their faces with the wireframe on top; strokes on the faces,
        # cut edges, a colour per face or the fields of the grid take its place (the wireframe
        # would swallow strokes and cut edges from afar, and hide the fields of the grid).
        view.show_wire = (number > 0 and not StageLines.marks_faces(data) and data.size_field is None
                          and data.grid_uv is None)
        for obj in context.selected_objects:
            obj.select_set(False)
        source.hide_set(True)
        view.select_set(True)
        context.view_layer.objects.active = view
        session.stage_lines.show(view, data)
        self.state.source = source
        self.state.view = view
        self.state.stage = identifier
        self.state.point_count = len(data.points)
        self.state.parts.clear()
        for budget, skipped in data.part_values:
            part = self.state.parts.add()
            part.budget = int(budget)
            part.skipped = bool(skipped)
        return view

    def refresh(self) -> None:
        """Computes the open stage again; errors go to the console (no operator to report them)."""
        session = AddonSession.current()
        if not self.is_open or self.state.stage not in _BY_ID or session is None or not session.core.ready:
            return
        try:
            self.show(self.state.stage, session)
        except CoreRunError as error:
            print(f"ArekoRemesher: {report_text(error.message)}")
        except CoreError as error:
            print(f"ArekoRemesher: {report_text('Internal error')}: {error}")

    def close(self) -> None:
        source = self.state.source
        session = AddonSession.current()
        if session is not None:
            session.stage_lines.clear()
        if self.state.view is not None:
            StageObject.remove(self.state.view)
        self.state.view = None
        self.state.stage = ""
        self.state.point_count = 0
        self.state.parts.clear()
        if source is not None:
            source.hide_set(False)
            for obj in self.context.selected_objects:
                obj.select_set(False)
            source.select_set(True)
            self.context.view_layer.objects.active = source
        self.state.source = None


class StageRefresh:
    """Computes an open stage object again after setting changes.

    Changes are collected and computed once, a moment after the last one, so
    dragging a slider does not compute the stage for every value. Only changes
    of settings that the shown stage reads count.
    """

    DELAY = 0.4  # seconds after the last change

    def __init__(self):
        self._scenes: set[str] = set()
        self._callback = self._run  # one bound method, so the timer can be found again

    def request(self, context, first_stage: int) -> None:
        view = InsightView(context)
        if not view.is_open or _BY_ID.get(view.state.stage, (-1, ""))[0] < first_stage:
            return
        self._scenes.add(context.scene.name)
        if bpy.app.timers.is_registered(self._callback):
            bpy.app.timers.unregister(self._callback)
        bpy.app.timers.register(self._callback, first_interval=self.DELAY)

    def cancel(self) -> None:
        if bpy.app.timers.is_registered(self._callback):
            bpy.app.timers.unregister(self._callback)
        self._scenes.clear()

    def _run(self) -> None:
        scenes, self._scenes = self._scenes, set()
        for window in bpy.context.window_manager.windows:
            if window.scene.name not in scenes:
                continue
            scenes.discard(window.scene.name)
            with bpy.context.temp_override(window=window):
                InsightView(bpy.context).refresh()
        return None


class AREKO_OT_show_stage(bpy.types.Operator):
    """Shows a stage as its own object"""

    bl_idname = "areko.show_stage"
    bl_label = "Show Stage"
    bl_options = {"REGISTER", "UNDO"}

    stage: EnumProperty(name="Stage", items=stage_items, translation_context=STAGE_CONTEXT)

    @classmethod
    def poll(cls, context):
        session = AddonSession.current()
        return session is not None and session.core.ready and InsightView(context).source() is not None

    def execute(self, context):
        try:
            InsightView(context).show(self.stage, AddonSession.current())
        except CoreRunError as error:
            self.report({"ERROR"}, report_text(error.message))
            return {"CANCELLED"}
        except CoreError as error:
            self.report({"ERROR"}, f"{report_text('Internal error')}: {error}")
            return {"CANCELLED"}
        return {"FINISHED"}


class AREKO_OT_close_stage(bpy.types.Operator):
    """Removes the stage object and shows the original again"""

    bl_idname = "areko.close_stage"
    bl_label = "Close Stage"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return InsightView(context).is_open

    def execute(self, context):
        InsightView(context).close()
        return {"FINISHED"}


CLASSES = (AREKO_PG_stage_part, AREKO_PG_insight, AREKO_OT_show_stage, AREKO_OT_close_stage)
