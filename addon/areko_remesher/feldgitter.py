# Grid of fields of a parameterization: the triangles of a stage mesh cut
# along the lines of the grid, into the fields of the grid cells and the lines
# between them. Fields and lines look the same on both sides of a cut; cells
# in two colours in turn would swap them wherever a cut carries a quarter turn.

from __future__ import annotations

import numpy as np

from .kern import StageData

# Colour numbers of the pieces: the grid lines, the fields between them, and
# the colour of triangles without grid coordinates.
LINE, FIELD, NO_FIELD = 0, 1, 2


class FieldGrid:
    """Cuts every triangle with grid coordinates into its pieces per field and grid line."""

    # Width of a grid line as a share of its cell; the field takes the rest. The lines
    # grow with the cells, so the grid reads the same from near and far.
    LINE_SHARE = 0.16
    # A triangle that covers more pieces than this (a degenerate one) stays whole.
    MAX_PIECES = 16384
    # Columns of a piece corner: position (0..2), normal (3..5), grid coordinates.
    _U, _V = 6, 7

    @classmethod
    def build(cls, data: StageData) -> StageData:
        """The stage with the cut mesh in place of its triangles.

        Every piece is a face of its own with its own points, a colour number as material and
        the normals of the stage mesh carried over, so that the grid is shaded like the mesh.
        """
        corners = data.corner_points.reshape(-1, 3)
        triangles = np.concatenate((data.points[corners].astype(np.float64), cls._point_normals(data)[corners],
                                    data.grid_uv.reshape(-1, 3, 2).astype(np.float64)), axis=2)
        uv = triangles[:, :, cls._U:]
        solved = np.isfinite(uv).all(axis=(1, 2))
        if data.face_parts is not None and data.grid_parts is not None:
            solved &= data.grid_parts[data.face_parts] != 0
        low, high = cls._reach(np.where(solved[:, None], uv.min(axis=1), 0.0),
                               np.where(solved[:, None], uv.max(axis=1), 0.0))
        covered = (high - low + 1.0).prod(axis=1)
        whole = ~solved | (covered == 1.0) | (covered > cls.MAX_PIECES)

        # Triangles in one strip along u and v (and those that stay whole) keep their shape.
        # A piece belongs to a field where both its strips do, and to a grid line otherwise.
        # A triangle in one strip takes the strip of its lowest coordinate, as the pieces of
        # the cut triangles do; one that covers too many pieces takes that of its centre.
        single = (covered[whole] == 1.0)[:, None]
        place = np.where(single, uv[whole].min(axis=1), uv[whole].mean(axis=1))
        place = np.where(solved[whole, None], place, 0.0)
        strips = cls._strip(np.clip(place, -1e15, 1e15)).astype(np.int64)
        colour = np.where(solved[whole], strips[:, 0] & strips[:, 1] & 1, NO_FIELD)
        pieces = [(triangles[whole], np.full(int(whole.sum()), 3, dtype=np.int64), colour)]
        # The others are cut into strips along u, the strips into pieces along v.
        three = np.full(int((~whole).sum()), 3, dtype=np.int64)
        for strip, strip_count, strip_u, _ in cls._split(triangles[~whole], three, cls._U, three):
            for piece, piece_count, strip_v, of_u in cls._split(strip, strip_count, cls._V, strip_u):
                pieces.append((piece, piece_count, of_u & strip_v & 1))

        counts = np.concatenate([count for _, count, _ in pieces])
        flat = np.concatenate([poly[np.arange(poly.shape[1])[None, :] < count[:, None]]
                               for poly, count, _ in pieces])
        normals = flat[:, 3:6]
        length = np.linalg.norm(normals, axis=1, keepdims=True)
        normals = np.divide(normals, length, out=np.zeros_like(normals), where=length > 0.0)
        none = np.empty(0, dtype=np.int32)
        return StageData(
            stage=data.stage,
            points=flat[:, :3].astype(np.float32),
            face_starts=np.concatenate(([0], np.cumsum(counts))).astype(np.int32),
            corner_points=np.arange(len(flat), dtype=np.int32),
            materials=np.concatenate([colour for _, _, colour in pieces]).astype(np.int32),
            face_groups=None, smoothing_groups=None, mirror_sides=None, face_parts=None,
            part_values=data.part_values,
            sharp_edges=np.empty((0, 2), dtype=np.int32),
            line_starts=np.zeros(1, dtype=np.int32), line_points=none, line_kinds=none,
            corner_normal=np.arange(len(flat), dtype=np.int32),
            normals=normals.astype(np.float32),
            colors=None, density=None, uv_sets=[],
        )

    @staticmethod
    def _point_normals(data: StageData) -> np.ndarray:
        """One normal per point of a triangle mesh: its triangle normals, weighted by area."""
        corners = data.corner_points.reshape(-1, 3)
        p = data.points[corners].astype(np.float64)
        face = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        normals = np.zeros((len(data.points), 3), dtype=np.float64)
        np.add.at(normals, corners.reshape(-1), np.repeat(face, 3, axis=0))
        return normals

    @classmethod
    def _strip(cls, values: np.ndarray) -> np.ndarray:
        """Number of the strip a grid coordinate lies in.

        Every grid coordinate is cut into strips in turn: strip 2n is the grid line n (as wide
        as LINE_SHARE, the line in its middle), strip 2n + 1 the fields between the lines n and
        n + 1. A strip begins at its lower bound and ends before the next one.
        """
        line = np.floor(values + 0.5)
        half = 0.5 * cls.LINE_SHARE
        return 2.0 * line + (values >= line + half) - (values < line - half)

    @classmethod
    def _bound(cls, strips: np.ndarray) -> np.ndarray:
        """Grid coordinate at which a strip begins."""
        line = np.floor(0.5 * strips)
        return line + np.where(strips == 2.0 * line, -0.5, 0.5) * cls.LINE_SHARE

    @classmethod
    def _reach(cls, lowest: np.ndarray, highest: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """First and last strip between two grid coordinates.

        A strip that the highest coordinate only touches at its beginning does not count.
        """
        first, last = cls._strip(lowest), cls._strip(highest)
        return first, np.maximum(last - (highest == cls._bound(last)), first)

    @classmethod
    def _split(cls, polygons: np.ndarray, counts: np.ndarray, column: int, carried: np.ndarray):
        """Cuts polygons into their parts per strip along one grid coordinate.

        `polygons` is (n, k, columns) with `counts` valid corners each; `carried` is one value
        per polygon that goes with its parts. Yields (parts, their corner counts, their strip
        numbers, their carried values), first the polygons that lie in one strip as they are,
        then the parts of the others strip by strip.
        """
        if len(polygons) == 0:
            return
        valid = np.arange(polygons.shape[1])[None, :] < counts[:, None]
        values = polygons[:, :, column]
        low, high = cls._reach(np.where(valid, values, np.inf).min(axis=1),
                               np.where(valid, values, -np.inf).max(axis=1))
        spans = (high - low).astype(np.int64) + 1
        single = spans == 1
        yield polygons[single], counts[single], low[single].astype(np.int64), carried[single]
        polygons, counts, carried, low, spans = (a[~single] for a in (polygons, counts, carried, low, spans))
        for offset in range(int(spans.max()) if len(spans) else 0):
            chosen = spans > offset
            strip = low[chosen] + offset
            part, count = cls._clip(polygons[chosen], counts[chosen], column, cls._bound(strip), True)
            part, count = cls._clip(part, count, column, cls._bound(strip + 1.0), False)
            kept = count >= 3
            yield part[kept], count[kept], strip[kept].astype(np.int64), carried[chosen][kept]

    @staticmethod
    def _clip(polygons: np.ndarray, counts: np.ndarray, column: int, bound: np.ndarray, above: bool):
        """Keeps of every polygon the side of the line `column` = `bound` (one bound per polygon).

        A corner on the line belongs to the kept side. Returns the polygons and their corner counts.
        """
        width = polygons.shape[1]
        index = np.arange(width)[None, :]
        valid = index < counts[:, None]
        following = np.take_along_axis(polygons, ((index + 1) % np.maximum(counts, 1)[:, None])[:, :, None], axis=1)
        here = polygons[:, :, column] - bound[:, None]
        there = following[:, :, column] - bound[:, None]
        if not above:
            here, there = -here, -there
        inside = here >= 0.0
        # An edge is cut where it passes the line between its two ends.
        crossing = valid & (inside != (there >= 0.0)) & (here != 0.0) & (there != 0.0)
        share = np.divide(here, here - there, out=np.zeros_like(here), where=crossing)
        cut = polygons + share[:, :, None] * (following - polygons)
        cut[:, :, column] = bound[:, None]
        both = np.empty((len(polygons), 2 * width, polygons.shape[2]), dtype=polygons.dtype)
        both[:, 0::2], both[:, 1::2] = polygons, cut
        keep = np.empty((len(polygons), 2 * width), dtype=bool)
        keep[:, 0::2], keep[:, 1::2] = valid & inside, crossing
        kept = keep.sum(axis=1)
        # Kept corners first, in their order.
        order = np.argsort(~keep, axis=1, kind="stable")[:, :max(int(kept.max()) if len(kept) else 0, 3)]
        return np.take_along_axis(both, order[:, :, None], axis=1), kept
