"""Build consume_plotfiles command lines for GRTeclyn wrapper episodes."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Literal, Sequence

from .config import REPO_ROOT, WRAPPER_ROOT
from .episode import Episode

ConsumerProfile = Literal["wormhole", "radial"]


def _strip_param_value(value: str) -> str:
    value = value.split("#", 1)[0].strip()
    if value.startswith('"') and value.endswith('"'):
        value = value[1:-1]
    return value


def _read_float_param(params_path: Path, key: str, default: float) -> float:
    try:
        for line in params_path.read_text(encoding="utf-8").splitlines():
            if "=" not in line:
                continue
            lhs, rhs = line.split("=", 1)
            if lhs.strip() == key:
                return float(_strip_param_value(rhs).split()[0])
    except (FileNotFoundError, ValueError, IndexError):
        pass
    return default


def _read_vector_param(params_path: Path, key: str, default: Sequence[float]) -> tuple[float, ...]:
    try:
        for line in params_path.read_text(encoding="utf-8").splitlines():
            if "=" not in line:
                continue
            lhs, rhs = line.split("=", 1)
            if lhs.strip() != key:
                continue
            values = tuple(float(item) for item in _strip_param_value(rhs).split())
            if values:
                return values
    except (FileNotFoundError, ValueError):
        pass
    return tuple(float(item) for item in default)


def resolve_consume_python() -> list[str]:
    """Return argv prefix for Python with visualization deps (yt)."""
    return [sys.executable]


def build_consume_command(
    episode: Episode,
    *,
    profile: ConsumerProfile = "radial",
    radii: Sequence[float] = (8.0, 16.0),
    n_points: int = 64,
    delete: bool = True,
    keep_last: int = 1,
    watch: bool = True,
    jobs: int = 4,
    frames: bool = True,
    keep_existing_frames: bool = False,
    stable_seconds: float | None = None,
    ftl_timeseries: bool = False,
    ftl_L: float | None = None,
) -> list[str]:
    """Return argv for streaming plotfile extraction into episode/small_data."""
    center = _read_vector_param(episode.params_path, "center", (0.0, 0.0, 0.0))
    if len(center) < 3:
        center = (*center, *([0.0] * (3 - len(center))))
    center = center[:3]
    l_full = _read_float_param(episode.params_path, "L_full", 40.0)
    # Frame/projection rendering is the dominant per-plotfile cost (12 yt renders
    # on a refined AMR grid, 10-70s each at max_level=3) and is *not* needed for
    # QD scoring -- only the FTL time-series + scalar metrics are.  Per-eval movies
    # are an inspection tool for promoted candidates, so the QD loop sets
    # GRTECLYN_FRAMES=0 to keep the consumer fast and stop the post-processing
    # backlog from starving the search.  Promote/inspection paths leave it unset.
    if os.environ.get("GRTECLYN_FRAMES", "").strip().lower() in {"0", "off", "no", "none", "false"}:
        frames = False
    zoom_env = os.environ.get("GRTECLYN_FRAMES_ZOOM", "").strip().lower()
    frame_zoom: float | None
    if zoom_env in {"", "none", "off", "full", "no", "0"}:
        frame_zoom = None
    else:
        frame_zoom = float(zoom_env)

    command = [
        *resolve_consume_python(),
        "-m",
        "grteclyn_wrapper.visualisation.process_wave.consume_plotfiles",
        "--data",
        str(episode.path),
        "--out",
        str(episode.small_data_dir),
        "--radii",
        *[str(r) for r in radii],
        "--n-points",
        str(n_points),
        "--center",
        *[f"{value:g}" for value in center],
        "--areal-radius",
        "-j",
        str(jobs),
        "--verbose",
    ]

    if watch:
        command.append("--watch")
    if delete:
        command.extend(["--delete", "--keep-last", str(keep_last)])
    if keep_existing_frames:
        command.append("--keep-existing-frames")
    if stable_seconds is not None:
        command.extend(["--stable-seconds", f"{float(stable_seconds):g}"])
    if ftl_timeseries:
        command.append("--ftl-timeseries")
        if ftl_L is not None:
            command.extend(["--ftl-l", f"{float(ftl_L):g}"])

    if profile == "wormhole":
        command.extend(
            [
                "--psi4",
                "--frames-fields",
                "K",
                "Weyl4_Re",
                "--frames-axis",
                "z",
                "--frames-corner",
                "--frames-out",
                str(episode.frames_dir),
                "--embedding",
                "--embedding-rmax",
                "5.0",
            ]
        )
    else:
        command.extend(
            [
                "--no-psi4",
                "--shell-fields",
                "chi",
                "lapse",
                "K",
            ]
        )
        if frames:
            # Which fields to render as slice-frame movies. Trace/gauge fields
            # (chi/lapse/K) are near-trivial for weak, momentum-carrying scalar
            # matter and look identical across candidates; the cloud and its
            # momentum live in phi/Pi, the shift (frame dragging), the
            # off-diagonal metric (shear) and rho_req. Override the default set
            # via $GRTECLYN_FRAMES_FIELDS (space-separated).
            env_fields = os.environ.get("GRTECLYN_FRAMES_FIELDS", "").split()
            frame_fields = env_fields if env_fields else ["chi", "lapse", "K"]
            projection_fields = os.environ.get("GRTECLYN_PROJECTION_FIELDS", "").split()
            projection_axes = os.environ.get("GRTECLYN_PROJECTION_AXES", "").split()
            command.extend(
                [
                    "--frames-fields",
                    *frame_fields,
                    "--frames-axis",
                    "z",
                    "--frames-center",
                    *[f"{value:g}" for value in center],
                    "--frames-out",
                    str(episode.frames_dir),
                ]
            )
            if frame_zoom is not None:
                command.extend(["--frames-zoom", f"{frame_zoom:g}"])
            if projection_fields and projection_axes:
                command.extend(
                    [
                        "--projection-fields",
                        *projection_fields,
                        "--projection-axes",
                        *projection_axes,
                        "--projection-method",
                        os.environ.get("GRTECLYN_PROJECTION_METHOD", "mip"),
                    ]
                )

    return command


def default_radii_for_example(example_name: str) -> tuple[float, ...]:
    if example_name == "RadialRecipe":
        return (4.0, 8.0)
    return (12.0, 16.0, 20.0, 24.0)
