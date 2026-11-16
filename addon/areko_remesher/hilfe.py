# Help behind the "?" beside every setting and button of the panels: one operator
# whose tooltip is the help text of its topic. English here; German in texte.py.

from __future__ import annotations

import bpy
from bpy.props import StringProperty

# Help text per topic (a setting name or a button), 2-4 sentences: what it is for and when to use it.
HELP = {
    # create quads
    "size": "Sets how fine the result gets. Count aims at about this many quads, Share takes a percentage of the "
            "input's faces, Edge Length gives the length of the quad edges in object units. Use Count for a fixed "
            "budget, e.g. a game asset, and Edge Length when quads should have the same size on several objects.",
    "exact_count": "Off, the number of quads may change a little to suit the shape. On, the count is kept strictly. "
                   "Turn it on when a fixed budget matters more than the flow, e.g. for a level of detail with a "
                   "strict limit.",
    "create_quads": "Computes a new quad mesh from the selected mesh object, with its modifiers applied; the "
                    "original mesh is not changed. Blender stays usable meanwhile, Esc stops the run. The result "
                    "is a new object named after the input. Shortcut in the 3D viewport: Ctrl Shift Q.",
    "cancel_quads": "Stops the running computation at once. Nothing is created and the input stays as it was. Use "
                    "it when a run takes too long or you want to change a setting first.",
    # density
    "curvature": "Makes the quads smaller where the surface bends strongly and larger on flat areas. Raise it for "
                 "figures and organic shapes, so that ears, fingers or folds get enough quads. Set 0 for quads of "
                 "about the same size everywhere, e.g. for a cloth simulation.",
    "attribute_density": "Reads the float point attribute \"density\" of the input: values up to 4 give smaller "
                         "quads, values down to 0.25 larger ones, 1 changes nothing. Use it when you build the "
                         "density yourself, e.g. with Geometry Nodes. By hand, Paint Density is easier.",
    "color_density": "Reads the first colour attribute of the input as a density map: red gives smaller quads, "
                     "cyan larger ones, white changes nothing. Turn it on to give single areas more or fewer quads, "
                     "e.g. the face of a figure. Paint Density turns it on by itself.",
    "density_paint": "Colour of the brush for Paint Density: 4 (red) gives the smallest quads, 0.25 (cyan) the "
                     "largest, 1 (white) changes nothing. You can change it while painting; the brush takes the new "
                     "colour at once.",
    "paint_density": "Switches the active object to Vertex Paint with the colour of Brush Density and turns on "
                     "Density from Color. A mesh without colours first gets a white colour attribute; otherwise its "
                     "first colour attribute is painted on, and its colours there are overwritten. Paint red where "
                     "you want more detail and cyan where fewer quads will do.",
    # fine control
    "points_per_circle": "How many quads go around a full round bend, e.g. around a cylinder. 0 takes the value "
                         "from Curvature Densifies. Raise it when round parts such as pipes or fingers look too "
                         "angular.",
    "density_ratio": "Limits how much larger the largest quads may be than the smallest. 0 takes the value from "
                     "Curvature Densifies, 1 gives all quads the same size. Raise it when large flat areas should "
                     "stay coarse while small details get fine quads.",
    "max_quad_ratio": "How long a quad may get compared to its width. Lower values keep the quads close to square, "
                      "higher values allow longer quads where the shape asks for it. The default suits most "
                      "models.",
    # edge flow
    "detect_creases": "Finds hard edges by the angle between neighbouring faces (Crease Angle) and lets edge loops "
                      "follow them. Keep it on for hard-surface and CAD models. Turn it off for smooth organic "
                      "shapes and scans, where small kinks should not pull loops.",
    "crease_angle": "With Detect Creases, faces that meet at a sharper angle than this count as a crease. Lower it "
                    "so that softer edges count too, e.g. on bevelled parts; raise it when loops land on small "
                    "kinks that do not matter.",
    "follow_materials": "Edge loops run along the borders between materials, so each material gets clean edges of "
                        "its own. Turn it on when the materials mark real parts, e.g. the panels of a hard-surface "
                        "model. It only acts when the object has at least two different materials.",
    "follow_face_groups": "Edge loops run along the borders between face sets, which you draw in Sculpt Mode. Turn "
                          "it on to place loops exactly where you want them, e.g. around an eye or along a panel "
                          "line, without touching the materials.",
    "follow_smoothing_groups": "Edge loops run along the borders between smoothing groups. Blender itself has none; "
                               "this reads the face attribute \"smoothing_groups\" that the results of this add-on "
                               "carry. Turn it on when you remesh such a result again and want to keep its "
                               "borders.",
    "follow_sharp_edges": "Edge loops run along the edges marked sharp in Edit Mode. It is the most direct way to "
                          "force a loop along a hard edge, e.g. on hard-surface and CAD models. Leave it off when the "
                          "sharp marks only serve the shading.",
    "follow_normal_splits": "Edge loops run along the edges where the normals split, e.g. custom normals of "
                            "imported game or CAD models. Use it when the hard edges exist only in the normals. On "
                            "a flat shaded model every edge splits, which usually gives a poor result.",
    "denoise": "Smooths small noise of the surface before the quads are computed; sharp edges stay. Keep it on for "
               "scans and sculpts, where it gives a calmer flow. Turn it off when fine details of a clean model "
               "should steer the flow.",
    "face_groups_to_materials": "Gives every face set of the selected object a new material in a colour of its own; "
                                "each face then uses the material of its face set. This shows the face sets "
                                "outside Sculpt Mode and lets Follow Material Borders use them. It needs at least "
                                "two face sets.",
    # guide curves
    "use_guide_curves": "Switches all guide curves of the list on or off at once. Turn it off to compare the result "
                        "with and without curves without losing the list.",
    "guide_curve_list": "Curve objects that steer the edge flow. Select curve objects and press + to add them, - "
                        "removes the chosen entry, and the checkbox of an entry switches that curve on or off. Draw "
                        "the curves on or close to the surface, e.g. around a neck or along a panel.",
    "curve": "The curve object of the chosen entry. Pick another one here to swap the curve and keep kind and "
             "width.",
    "kind": "How the curve steers the result: Band makes the quads beside it follow its direction, Edge Loop puts "
            "an edge loop along it, Crease does the same and keeps that loop sharp, Off ignores it. Use Edge Loop "
            "for exact loops, e.g. around a mouth, and Band for a general direction, e.g. along a limb.",
    "automatic_width": "On, the width of the curve's influence follows from its length. Turn it off to set Width "
                       "yourself, e.g. when a Band should steer a wider or narrower stripe.",
    "width": "How far beside the curve the quads follow its direction, in object units. It only counts when "
             "Automatic Width is off. Widen it for a broad, calm flow; narrow it when the curve should act only "
             "close to itself.",
    "new_kind": "The kind that curves get when you add them with +. Set it before you add several curves of the "
                "same kind; each entry can still be changed afterwards.",
    # borders
    "follow_borders": "Edge loops run along the open borders of the mesh, e.g. sleeve openings, hems or the rim of "
                      "a hole. Keep it on for clothing and other open shapes; turn it off only when the shape of "
                      "the border should not steer the flow.",
    "freeze_borders": "The points of the open borders stay exactly where they are, and the size of the quads then "
                      "follows the spacing of these points, not the count. Use it when the result must fit another "
                      "mesh along its border, e.g. a clothing piece that joins a body.",
    "border_density_factor": "Only with Freeze Borders: scales the size of the quads that comes from the spacing of "
                             "the frozen border points. 1 keeps that spacing, values above 1 give larger quads, "
                             "values below 1 smaller ones.",
    # seams
    "weld_seams": "Joins open seams whose border points lie exactly on top of each other, so the mesh is closed "
                  "there. Turn it on for models that look closed but are split along seams, e.g. some imports; "
                  "otherwise each seam counts as an open border.",
    "weld_as_crease": "Seams joined by Weld Open Seams count as creases: edge loops follow them and they stay "
                      "sharp. Use it when the seams are real edges, e.g. where the panels of a CAD model meet; "
                      "leave it off when they were only a technical split.",
    "merge_cad_groups": "Merges narrow strips between face groups into their neighbours, such as the small fillets "
                        "of CAD exports, so they do not pull edge loops of their own. It acts on the borders that "
                        "Edge Flow follows (materials, face groups, smoothing groups).",
    # symmetry
    "mirror_axes": "Computes one half and mirrors it across the chosen planes, so the result is exactly symmetric. "
                   "X suits a figure that faces forward; several axes can be combined. Mirror at Object Center "
                   "decides where the planes lie.",
    "mirror_at_center": "On, the mirror plane runs through the middle of the object's bounding box; off, through "
                        "its origin. Turn it off when the object is symmetric about its origin but its outline is "
                        "not, e.g. a figure holding something in one hand.",
    "detect_symmetry": "Finds the mirror plane of every part by itself and computes one half of each, without an "
                       "axis. Use it for parts that are symmetric but not aligned with the axes, or for several "
                       "parts with planes of their own. It compares faces on both sides, so it works best on "
                       "modelled meshes.",
    "freeze_half": "For a mesh that is only one half, e.g. for a Mirror modifier: the open border on the mirror "
                   "plane stays in place, so both halves still fit together, and the result stays a half. Choose "
                   "the axis under Mirror on; it only makes sense together with Freeze Borders.",
    # output
    "hide_input": "Hides the input object as soon as the result is there, so the two never lie on top of each "
                  "other. Turn it off to keep both visible; then move one aside to compare them.",
    "split_at_creases": "Splits the points of the result along its creases, so the faces on either side no longer "
                        "share them. Use it for exports where hard edges must be real splits, e.g. for baking or a "
                        "game engine; leave it off when you keep modelling the result.",
    "creases_to_edge_crease": "Gives the creases of the result full edge crease, so a Subdivision Surface modifier "
                              "keeps them sharp. Turn it on for hard-surface models that you subdivide afterwards.",
    "clean_output": "Tidies the result: it dissolves stray edges and separates points where the mesh is not "
                    "manifold. Turn it on when the result goes to tools that need a clean mesh, e.g. a simulation "
                    "or 3D printing.",
    "keep_colors": "Carries the colour attribute of the input over to the result. Turn it on when the colours "
                   "matter later, e.g. painted vertex colours; leave it off when they only served as a density map.",
    "thread_count": "How many processor cores compute the separate pieces of the mesh side by side; 0 uses all. "
                    "It only speeds up meshes made of several pieces, and the result is the same for every number. "
                    "Lower it when other programs should stay responsive.",
    "remember_preparation": "Keeps the prepared mesh after a run. If you then change only settings that act later, "
                            "e.g. Curvature Densifies or the output, the next run skips the preparation and is "
                            "faster; the result is the same. Turn it off to free the memory with very large "
                            "meshes.",
    "reset": "Sets all settings of the panel back to their defaults, including Use Guide Curves and Kind for New "
             "Curves; the guide curve list stays. Use it for a clean start when a mix of settings gives odd "
             "results.",
    # insight
    "show_stage": "Computes up to the chosen stage and shows it as an object of its own in place of the input, "
                  "which is hidden meanwhile. Use it to see why a result looks as it does, e.g. Feature Lines shows "
                  "the creases found and Size Field where the quads get small. Changed settings update the shown "
                  "stage after a moment.",
    "close_stage": "Removes the stage object and shows the input again. Create Quads does this by itself before it "
                   "starts.",
}


class AREKO_OT_help(bpy.types.Operator):
    """Help"""

    bl_idname = "areko.help"
    bl_label = "Help"
    bl_options = {"INTERNAL"}

    topic: StringProperty(options={"HIDDEN", "SKIP_SAVE"})

    @classmethod
    def description(cls, context, properties):
        text = HELP.get(properties.topic)
        return bpy.app.translations.pgettext_tip(text) if text else ""

    def execute(self, context):
        # The button only carries its tooltip: a click changes nothing and leaves no undo step.
        return {"CANCELLED"}


CLASSES = (AREKO_OT_help,)
