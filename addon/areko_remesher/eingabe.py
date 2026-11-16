# Input of a run from a Blender object: the evaluated mesh (modifiers applied)
# in object coordinates, in Blender's point, face and corner order.

from __future__ import annotations

import bpy
import numpy as np

from .kern import MeshInput

# Point attribute read as density factor.
DENSITY_ATTRIBUTE = "density"
FACE_SET_ATTRIBUTE = ".sculpt_face_set"
# Face attribute of smoothing-group bit masks (written by ausgabe.ResultObject).
SMOOTHING_GROUPS_ATTRIBUTE = "smoothing_groups"


class BlenderMeshInput:
    """Builds the core input from one mesh object."""

    @classmethod
    def from_object(cls, obj: bpy.types.Object, depsgraph) -> MeshInput:
        evaluated = obj.evaluated_get(depsgraph)
        mesh = evaluated.to_mesh(preserve_all_data_layers=True, depsgraph=depsgraph)
        try:
            data = cls.from_mesh(mesh, [slot.material for slot in obj.material_slots])
        finally:
            evaluated.to_mesh_clear()
        # A mesh that keeps its shape keys keeps one normal per corner.
        data.share_normals = not cls._keeps_shape_keys(obj)
        return data

    @staticmethod
    def _keeps_shape_keys(obj: bpy.types.Object) -> bool:
        """Shape keys beside the basis, and no active modifier that would replace the mesh."""
        keys = getattr(obj.data, "shape_keys", None)
        if keys is None or len(keys.key_blocks) < 2:
            return False
        return not any(m.show_viewport or m.show_render for m in obj.modifiers if m.type != "ARMATURE")

    @staticmethod
    def _read(collection, prop: str, count: int, dtype, width: int = 1) -> np.ndarray:
        values = np.empty(count * width, dtype=dtype)
        collection.foreach_get(prop, values)
        return values.reshape(-1, width) if width > 1 else values

    @classmethod
    def from_mesh(cls, mesh: bpy.types.Mesh, materials: list) -> MeshInput:
        n_points, n_faces, n_corners = len(mesh.vertices), len(mesh.polygons), len(mesh.loops)
        corner_points = cls._read(mesh.loops, "vertex_index", n_corners, np.int32)
        data = MeshInput(
            points=cls._read(mesh.vertices, "co", n_points, np.float32, 3),
            face_sizes=cls._read(mesh.polygons, "loop_total", n_faces, np.int32),
            corner_points=corner_points,
        )
        attributes = mesh.attributes

        # Normals of every corner, also on smooth meshes.
        data.corner_normals = cls._read(mesh.corner_normals, "vector", n_corners, np.float32, 3)

        if len(mesh.color_attributes):
            first = mesh.color_attributes[0]
            per_point = first.domain == "POINT"
            colors = cls._read(first.data, "color_srgb", n_points if per_point else n_corners, np.float32, 4)
            data.corner_colors = colors[corner_points] if per_point else colors

        distinct = []
        for material in materials:
            if material is not None and all(material is not m for m in distinct):
                distinct.append(material)
        index = attributes.get("material_index")
        if len(distinct) > 1 and index is not None:
            data.materials = cls._read(index.data, "value", n_faces, np.int32)

        face_sets = attributes.get(FACE_SET_ATTRIBUTE)
        if face_sets is not None and face_sets.domain == "FACE" and face_sets.data_type == "INT":
            data.face_groups = cls._read(face_sets.data, "value", n_faces, np.int32)

        # Smoothing groups only when the mesh carries them (e.g. a result of this add-on); none are derived
        # from sharp edges, as the reference add-on's FBX export ("Normals Only") never writes any.
        smoothing = attributes.get(SMOOTHING_GROUPS_ATTRIBUTE)
        if smoothing is not None and smoothing.domain == "FACE" and smoothing.data_type == "INT":
            data.smoothing_groups = cls._read(smoothing.data, "value", n_faces, np.int32)

        sharp = attributes.get("sharp_edge")
        if sharp is not None and sharp.domain == "EDGE":
            flags = cls._read(sharp.data, "value", len(mesh.edges), bool)
            if flags.any():
                edges = cls._read(mesh.edges, "vertices", len(mesh.edges), np.int32, 2)
                data.sharp_edges = edges[flags]

        density = attributes.get(DENSITY_ATTRIBUTE)
        if density is not None and density.domain == "POINT" and density.data_type == "FLOAT":
            data.point_density = cls._read(density.data, "value", n_points, np.float32)

        for layer in mesh.uv_layers:
            data.uv_sets.append((layer.name, cls._read(layer.uv, "vector", n_corners, np.float32, 2)))
        return data
