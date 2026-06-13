from __future__ import annotations

from pathlib import Path
from typing import Sequence


def _cleanup_existing_frames(frames_out_dir: str, fields: Sequence[str], axis: str, verbose: bool) -> None:
    """
    Remove existing PNG frames (and movie mp4) for the requested field+axis outputs.

    This mirrors the behavior of visualize/__main__.py which clears frames before
    creating a new animation run, but is scoped to the fields requested via
    --frames-fields (so it won't touch other outputs).
    """
    base = Path(frames_out_dir)
    for fld in fields:
        out_dir = base / f"{fld}_{axis}"
        frames_dir = out_dir / "frames"
        if frames_dir.is_dir():
            for p in frames_dir.glob("*.png"):
                try:
                    p.unlink()
                except FileNotFoundError:
                    pass
            if verbose:
                print(f"[clean] cleared frames in {frames_dir}")
        # Also remove a previously stitched movie for this field/axis if present.
        movie = out_dir / f"movie_{fld}_{axis}.mp4"
        if movie.exists():
            try:
                movie.unlink()
                if verbose:
                    print(f"[clean] removed {movie}")
            except FileNotFoundError:
                pass


def _cleanup_projection_frames(frames_out_dir: str, fields: Sequence[str], axes: Sequence[str], verbose: bool) -> None:
    """Remove existing projection PNG frames for requested field+axis outputs."""
    base = Path(frames_out_dir)
    for fld in fields:
        for axis in axes:
            out_dir = base / f"{fld}_proj_{axis}"
            frames_dir = out_dir / "frames"
            if frames_dir.is_dir():
                for p in frames_dir.glob("*.png"):
                    try:
                        p.unlink()
                    except FileNotFoundError:
                        pass
                if verbose:
                    print(f"[clean] cleared projection frames in {frames_dir}")
            movie = out_dir / f"movie_{fld}_proj_{axis}.mp4"
            if movie.exists():
                try:
                    movie.unlink()
                    if verbose:
                        print(f"[clean] removed {movie}")
                except FileNotFoundError:
                    pass


def _cleanup_embedding_frames(frames_out_dir: str, verbose: bool) -> None:
    """Remove existing embedding PNG frames and movie."""
    base = Path(frames_out_dir) / "embedding" / "frames"
    if base.is_dir():
        for p in base.glob("*.png"):
            try:
                p.unlink()
            except FileNotFoundError:
                pass
        if verbose:
            print(f"[clean] cleared embedding frames in {base}")
    movie = Path(frames_out_dir) / "embedding" / "movie_embedding.mp4"
    if movie.exists():
        try:
            movie.unlink()
        except FileNotFoundError:
            pass
