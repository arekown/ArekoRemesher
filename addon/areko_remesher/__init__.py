# ArekoRemesher - Blender extension: registration, sidebar panels and shortcut.

from __future__ import annotations

import bpy

from . import bedienung, einblick, hilfe, kurven
from .ausgabe import StageLines
from .einstellungen import SIZE_FIELD, settings_group
from .sitzung import AddonSession
from .texte import PANEL_CONTEXT, STAGE_CONTEXT, Translations, decimal, grouped, iface

TAB = "Areko"
# Shortcut of "Create Quads" in the 3D viewport (free in Blender's default keymap there).
SHORTCUT = {"type": "Q", "ctrl": True, "shift": True}


class _SidebarPanel:
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = TAB
    bl_translation_context = PANEL_CONTEXT


def _help_button(layout, topic: str) -> None:
    """The "?" whose tooltip explains `topic` (see hilfe.HELP)."""
    layout.operator(hilfe.AREKO_OT_help.bl_idname, text="", icon="QUESTION", emboss=False).topic = topic


def _helped(layout, topic: str):
    """One row: the element, drawn into the returned layout, and its "?" on the right."""
    row = layout.row(align=True)
    element = row.row(align=True)
    _help_button(row, topic)
    return element


def _hint(layout, *lines: str) -> None:
    """Greyed lines of short help below a setting."""
    column = layout.column(align=True)
    column.active = False
    for line in lines:
        column.label(text=line)


# Logo, name, core version and core status (comment, not a docstring: no tooltip).
class AREKO_PT_header(_SidebarPanel, bpy.types.Panel):
    bl_idname = "AREKO_PT_header"
    bl_label = "ArekoRemesher"
    bl_options = {"HIDE_HEADER"}
    bl_order = 0

    def draw(self, context):
        session = AddonSession.current()
        layout = self.layout
        row = layout.row(align=False)
        row.template_icon(icon_value=session.logo_icon, scale=3.0)
        text = row.column(align=True)
        text.label(text="ArekoRemesher", translate=False)
        if session.core.ready:
            text.label(text=f"{iface('Version')} {session.core.version}", translate=False)
            layout.label(text="Core ready", icon="CHECKMARK")
        else:
            layout.label(text="Core not available", icon="ERROR")


# Main action: size of the quads, the button, progress while computing and the last
# result (comment, not a docstring: no tooltip).
class AREKO_PT_create_quads(_SidebarPanel, bpy.types.Panel):
    bl_idname = "AREKO_PT_create_quads"
    bl_label = "Create Quads"
    bl_order = 1

    def draw(self, context):
        settings = context.scene.areko_remesher
        session = AddonSession.current()
        layout = self.layout
        size = _helped(layout, "size").split(factor=0.5, align=True)
        size.prop(settings, "size_kind", text="")
        size.prop(settings, SIZE_FIELD[settings.size_kind], text="")
        _helped(layout, "exact_count").prop(settings, "exact_count")
        layout.separator(type="LINE")
        if bedienung.flat_with_normal_splits(context, context.active_object):
            warning = layout.column(align=True)
            warning.label(text="Flat shading and", icon="ERROR")
            warning.label(text="Follow Normal Splits:", icon="BLANK1")
            warning.label(text="usually a poor result", icon="BLANK1")
        run = session.run
        if run is not None:
            layout.progress(factor=run.job.fraction, type="BAR", text=run.progress_text())
            _helped(layout, "cancel_quads").operator(bedienung.AREKO_OT_cancel_quads.bl_idname, icon="CANCEL")
            return
        big = layout.row()
        big.scale_y = 1.6
        row = _helped(big, "create_quads")
        row.enabled = session.core.ready
        row.operator(bedienung.AREKO_OT_create_quads.bl_idname, depress=True)
        summary = session.last_result
        if summary is not None and summary.name in bpy.data.objects:
            box = layout.column(align=True)
            box.active = False
            box.label(text=iface("{} Quads · {} Triangles").format(grouped(summary.quads),
                                                                    grouped(summary.triangles)), translate=False)
            box.label(text=iface("Time: {} s").format(decimal(summary.seconds, 1)), translate=False)


# Density: curvature, density maps and the density brush (comment, not a docstring: no tooltip).
class AREKO_PT_density(_SidebarPanel, bpy.types.Panel):
    bl_idname = "AREKO_PT_density"
    bl_label = "Density"
    bl_order = 2

    def draw(self, context):
        settings = context.scene.areko_remesher
        layout = self.layout
        _helped(layout, "curvature").prop(settings, "curvature", slider=True)
        _hint(layout, "0 % = the same size everywhere")
        _helped(layout, "attribute_density").prop(settings, "attribute_density")
        _helped(layout, "color_density").prop(settings, "color_density")
        layout.separator(type="LINE")
        _helped(layout, "density_paint").prop(settings, "density_paint", slider=True)
        _hint(layout, "4 = finer · 0.25 = coarser")
        _helped(layout, "paint_density").operator(bedienung.AREKO_OT_paint_density.bl_idname, icon="BRUSH_DATA")


# Fine control of the size field (comment, not a docstring: no tooltip).
class AREKO_PT_density_fine(_SidebarPanel, bpy.types.Panel):
    bl_idname = "AREKO_PT_density_fine"
    bl_label = "Fine Control"
    bl_parent_id = "AREKO_PT_density"
    bl_options = {"DEFAULT_CLOSED"}

    def draw(self, context):
        settings = context.scene.areko_remesher
        layout = self.layout
        for name in ("points_per_circle", "density_ratio", "max_quad_ratio"):
            _helped(layout, name).prop(settings, name)


# Edge flow: how the result follows the surface (comment, not a docstring: no tooltip).
class AREKO_PT_edge_flow(_SidebarPanel, bpy.types.Panel):
    bl_idname = "AREKO_PT_edge_flow"
    bl_label = "Edge Flow"
    bl_order = 3

    def draw(self, context):
        settings = context.scene.areko_remesher
        layout = self.layout
        for name in ("detect_creases", "crease_angle", "follow_materials", "follow_face_groups",
                     "follow_smoothing_groups", "follow_sharp_edges", "follow_normal_splits", "denoise"):
            _helped(layout, name).prop(settings, name)
        _helped(layout, "face_groups_to_materials").operator(bedienung.AREKO_OT_face_groups_to_materials.bl_idname,
                                                             icon="MATERIAL")


# Guide curves under the edge flow: list of curve objects with kind, on/off and
# width (comment, not a docstring: no tooltip).
class AREKO_PT_guide_curves(_SidebarPanel, bpy.types.Panel):
    bl_idname = "AREKO_PT_guide_curves"
    bl_label = "Guide Curves"
    bl_parent_id = "AREKO_PT_edge_flow"
    bl_order = 0

    def draw_header(self, context):
        self.layout.prop(context.scene.areko_curves, "use", text="")

    def draw_header_preset(self, context):
        _help_button(self.layout, "use_guide_curves")

    def draw(self, context):
        curves = context.scene.areko_curves
        layout = self.layout
        layout.active = curves.use
        row = layout.row()
        row.template_list("AREKO_UL_guide_curves", "", curves, "entries", curves, "index", rows=3)
        buttons = row.column(align=True)
        buttons.operator(kurven.AREKO_OT_add_guide_curves.bl_idname, icon="ADD", text="")
        buttons.operator(kurven.AREKO_OT_remove_guide_curve.bl_idname, icon="REMOVE", text="")
        buttons.separator()
        _help_button(buttons, "guide_curve_list")
        if 0 <= curves.index < len(curves.entries):
            entry = curves.entries[curves.index]
            _helped(layout, "curve").prop(entry, "curve", text="")
            _helped(layout, "kind").prop(entry, "kind", text="")
            _helped(layout, "automatic_width").prop(entry, "automatic_width")
            width = _helped(layout, "width")
            width.active = not entry.automatic_width
            width.prop(entry, "width")
        layout.label(text="Kind for New Curves")
        _helped(layout, "new_kind").prop(curves, "new_kind", text="")


# Open borders under the edge flow (comment, not a docstring: no tooltip).
class AREKO_PT_borders(_SidebarPanel, bpy.types.Panel):
    bl_idname = "AREKO_PT_borders"
    bl_label = "Borders"
    bl_parent_id = "AREKO_PT_edge_flow"
    bl_options = {"DEFAULT_CLOSED"}
    bl_order = 1

    def draw(self, context):
        settings = context.scene.areko_remesher
        layout = self.layout
        for name in ("follow_borders", "freeze_borders", "border_density_factor"):
            _helped(layout, name).prop(settings, name)


# Seams and CAD meshes under the edge flow (comment, not a docstring: no tooltip).
class AREKO_PT_seams(_SidebarPanel, bpy.types.Panel):
    bl_idname = "AREKO_PT_seams"
    bl_label = "Seams"
    bl_parent_id = "AREKO_PT_edge_flow"
    bl_options = {"DEFAULT_CLOSED"}
    bl_order = 2

    def draw(self, context):
        settings = context.scene.areko_remesher
        layout = self.layout
        for name in ("weld_seams", "weld_as_crease", "merge_cad_groups"):
            _helped(layout, name).prop(settings, name)


# Symmetry: mirror axes and how the plane is found (comment, not a docstring: no tooltip).
class AREKO_PT_symmetry(_SidebarPanel, bpy.types.Panel):
    bl_idname = "AREKO_PT_symmetry"
    bl_label = "Symmetry"
    bl_order = 4

    def draw(self, context):
        settings = context.scene.areko_remesher
        layout = self.layout
        split = _helped(layout, "mirror_axes").split(factor=0.45)
        split.label(text="Mirror on")
        axes = split.row(align=True)
        for axis in ("x", "y", "z"):
            axes.prop(settings, f"mirror_{axis}", text=axis.upper(), toggle=True, translate=False)
        for name in ("mirror_at_center", "detect_symmetry", "freeze_half"):
            _helped(layout, name).prop(settings, name)


# Output: what the result keeps, the cores, the reset and the shortcut (comment, not a docstring: no tooltip).
class AREKO_PT_output(_SidebarPanel, bpy.types.Panel):
    bl_idname = "AREKO_PT_output"
    bl_label = "Output"
    bl_order = 5

    def draw(self, context):
        settings = context.scene.areko_remesher
        layout = self.layout
        for name in ("hide_input", "split_at_creases", "creases_to_edge_crease", "clean_output", "keep_colors",
                     "thread_count", "remember_preparation"):
            _helped(layout, name).prop(settings, name)
        layout.separator(type="LINE")
        _helped(layout, "reset").operator(bedienung.AREKO_OT_reset_settings.bl_idname, icon="LOOP_BACK")
        if _keymap_items:
            _hint(layout, iface("Shortcut: {}").format(_keymap_items[0][1].to_string()))


# Stage insight: one stage as its own object. A comment, not a docstring, so the
# panel gets no tooltip.
class AREKO_PT_insight(_SidebarPanel, bpy.types.Panel):
    bl_idname = "AREKO_PT_insight"
    bl_label = "Insight"
    bl_options = {"DEFAULT_CLOSED"}
    bl_order = 6

    def draw(self, context):
        layout = self.layout
        hint = layout.row()
        hint.active = False
        hint.label(text="Shows a stage as its own object.")
        split = _helped(layout, "show_stage").split(factor=0.35)
        split.label(text="Stage")
        split.operator_menu_enum(einblick.AREKO_OT_show_stage.bl_idname, "stage",
                                 text=einblick.stage_label(context.scene.areko_insight.stage),
                                 text_ctxt=STAGE_CONTEXT)
        insight = context.scene.areko_insight
        if insight.view is not None:
            count = layout.row()
            count.active = False
            count.label(text=iface("{} Points").format(grouped(insight.point_count)), translate=False)
            self._draw_parts(layout, insight)
        _helped(layout, "close_stage").operator(einblick.AREKO_OT_close_stage.bl_idname)

    MAX_PARTS_SHOWN = 8

    @classmethod
    def _draw_parts(cls, layout, insight) -> None:
        if not len(insight.parts):
            return
        column = layout.column(align=True)
        column.active = False
        column.label(text=iface("{} Parts").format(grouped(len(insight.parts))), translate=False)
        for i, part in enumerate(insight.parts[:cls.MAX_PARTS_SHOWN]):
            if part.skipped:
                text = iface("Part {}: skipped").format(i + 1)
            else:
                text = iface("Part {}: {} Quads").format(i + 1, grouped(part.budget))
            column.label(text=text, translate=False)


_PANELS = (AREKO_PT_header, AREKO_PT_create_quads, AREKO_PT_density, AREKO_PT_density_fine, AREKO_PT_edge_flow,
           AREKO_PT_guide_curves, AREKO_PT_borders, AREKO_PT_seams, AREKO_PT_symmetry, AREKO_PT_output,
           AREKO_PT_insight)
_translations = Translations(__name__)
_keymap_items: list = []  # (keymap, item) of the shortcut


def _classes(session: AddonSession) -> tuple:
    return (session.settings_class, *bedienung.CLASSES, *kurven.CLASSES, *einblick.CLASSES, *hilfe.CLASSES, *_PANELS)


def _add_shortcut() -> None:
    keyconfig = bpy.context.window_manager.keyconfigs.addon
    if keyconfig is None:  # Blender without user interface
        return
    keymap = keyconfig.keymaps.new(name="3D View", space_type="VIEW_3D")
    item = keymap.keymap_items.new(bedienung.AREKO_OT_create_quads.bl_idname, value="PRESS", **SHORTCUT)
    _keymap_items.append((keymap, item))


def _remove_shortcut() -> None:
    for keymap, item in _keymap_items:
        keymap.keymap_items.remove(item)
    _keymap_items.clear()


def register():
    session = AddonSession.open()
    _translations.register()
    session.settings_class = settings_group(session.library.default_settings() if session.library else None)
    session.stage_refresh = einblick.StageRefresh()
    session.stage_lines = StageLines()
    for cls in _classes(session):
        bpy.utils.register_class(cls)
    bpy.types.Scene.areko_remesher = bpy.props.PointerProperty(type=session.settings_class)
    bpy.types.Scene.areko_insight = bpy.props.PointerProperty(type=einblick.AREKO_PG_insight)
    kurven.register_properties()
    _add_shortcut()


def unregister():
    session = AddonSession.current()
    if session.run is not None:
        session.run.job.cancel()
    _remove_shortcut()
    session.stage_refresh.cancel()
    session.stage_lines.clear()
    kurven.unregister_properties()
    del bpy.types.Scene.areko_insight
    del bpy.types.Scene.areko_remesher
    for cls in reversed(_classes(session)):
        bpy.utils.unregister_class(cls)
    _translations.unregister()
    AddonSession.close()
