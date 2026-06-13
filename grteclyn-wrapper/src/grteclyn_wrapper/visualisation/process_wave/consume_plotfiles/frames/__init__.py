"""Slice, projection, and embedding frame rendering."""

from __future__ import annotations

from .cleanup import (
    _cleanup_embedding_frames,
    _cleanup_existing_frames,
    _cleanup_projection_frames,
)
from .embedding import _render_embedding_frame
from .projection import _render_projection_frame
from .slice import _render_slice_frame
from .zlim import _lock_frame_zlims_from_plotfile

__all__ = [
    "_cleanup_embedding_frames",
    "_cleanup_existing_frames",
    "_cleanup_projection_frames",
    "_lock_frame_zlims_from_plotfile",
    "_render_embedding_frame",
    "_render_projection_frame",
    "_render_slice_frame",
]
