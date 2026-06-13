from __future__ import annotations

from typing import Sequence

import numpy as np

from ..config import _FRAME_BUFF_CAP, _FRAME_SAMPLES_PER_CELL


def _resolve_frame_physics_center(
    ds,
    axis: str,
    coord: float | None,
    zoom: float | None,
    center_xyz: Sequence[float] | None,
    corner: bool,
) -> list[float]:
    """Match visualize/__main__.py center logic; auto throat corner when unset."""
    mid_x = float((ds.domain_right_edge[0] + ds.domain_left_edge[0]) / 2.0)
    mid_y = float((ds.domain_right_edge[1] + ds.domain_left_edge[1]) / 2.0)
    physics_center = [mid_x, mid_y, 0.0]

    if corner and zoom is not None:
        slice_plane_val = 0.0 if coord is None else float(coord)
        if center_xyz is not None:
            origin = np.array(center_xyz, dtype=float)
        else:
            origin = np.array(
                [mid_x - float(zoom) / 2.0, mid_y - float(zoom) / 2.0, slice_plane_val],
                dtype=float,
            )
        w = float(zoom)
        if axis == "z":
            physics_center = [origin[0] + w / 2.0, origin[1] + w / 2.0, slice_plane_val]
        elif axis == "y":
            physics_center = [origin[0] + w / 2.0, slice_plane_val, origin[2] + w / 2.0]
        elif axis == "x":
            physics_center = [slice_plane_val, origin[1] + w / 2.0, origin[2] + w / 2.0]
    elif center_xyz is not None:
        physics_center = [float(center_xyz[0]), float(center_xyz[1]), float(center_xyz[2])]

    if coord is not None:
        if axis == "z":
            physics_center[2] = float(coord)
        elif axis == "y":
            physics_center[1] = float(coord)
        elif axis == "x":
            physics_center[0] = float(coord)

    return physics_center


def _frame_buff_size(ds, zoom: float | None) -> int:
    """Pixels across the FRB: ``_FRAME_SAMPLES_PER_CELL`` samples per simulation
    cell in the zoom window (capped at ``_FRAME_BUFF_CAP``)."""
    if zoom is None:
        zoom = float(ds.domain_width[0])
    le = np.array(ds.domain_left_edge.d, dtype=float)
    re = np.array(ds.domain_right_edge.d, dtype=float)
    dims = np.array(ds.domain_dimensions, dtype=float) * (2 ** int(ds.index.max_level))
    dx = (re - le) / np.maximum(dims, 1.0)
    cells = max(int(round(float(zoom) / min(dx[0], dx[1]))), 64)
    return int(min(max(cells * _FRAME_SAMPLES_PER_CELL, 256), _FRAME_BUFF_CAP))

