"""Build consume_plotfiles command lines for GRTeclyn wrapper episodes."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Literal, Mapping, Sequence

from .config import REPO_ROOT, WRAPPER_ROOT
from .episode import Episode

ConsumerProfile = Literal["wormhole", "radial"]


def _env_flag(name: str) -> bool:
    return os.environ.get(name, "").strip().lower() in {"1", "on", "yes", "true"}


def _central_timeseries_enabled(explicit: bool | None) -> bool:
    if explicit is not None:
        return explicit
    return _env_flag("GRTECLYN_CENTRAL_TIMESERIES")


def _confinement_timeseries_enabled(explicit: bool | None, ftl_timeseries: bool) -> bool:
    """Matter-confinement stream: the trustworthy dispersal detector.

    Default-on for any matter FTL run (it is a cheap covering-grid moment) so we
    always get the rms_radius / confined_frac curves; ``GRTECLYN_CONFINEMENT=0``
    disables it, and an explicit argument overrides both.
    """
    if explicit is not None:
        return explicit
    if os.environ.get("GRTECLYN_CONFINEMENT", "").strip().lower() in {
        "0", "off", "no", "none", "false",
    }:
        return False
    return bool(ftl_timeseries) or _env_flag("GRTECLYN_CONFINEMENT")


def _central_ball_enabled() -> bool:
    return _env_flag("GRTECLYN_CENTRAL_BALL")


def _central_radial_enabled() -> bool:
    return _env_flag("GRTECLYN_CENTRAL_RADIAL")


def _incremental_score_enabled(explicit: bool) -> bool:
    if explicit:
        return True
    return _env_flag("GRTECLYN_INCREMENTAL_SCORE")


def _splash_early_term_enabled() -> bool:
    return _env_flag("GRTECLYN_SPLASH_EARLY_TERM")


def _psi4_enabled() -> bool:
    return _env_flag("GRTECLYN_PSI4")


def consumer_frames_enabled() -> bool:
    """Whether the plotfile consumer renders PNG frames and projection panels.

    Disabled when ``GRTECLYN_FRAMES=0`` or ``GRTECLYN_CONSUMER_DRAIN=1`` (fast
    backlog drain: Psi4 + FTL metrics only, ~10x faster on HQ AMR dumps).
    """
    if _env_flag("GRTECLYN_CONSUMER_DRAIN"):
        return False
    if os.environ.get("GRTECLYN_FRAMES", "").strip().lower() in {
        "0",
        "off",
        "no",
        "none",
        "false",
    }:
        return False
    return True


def post_run_frames_enabled(*, explicit: bool | None = None) -> bool:
    """Default frame rendering for post-GPU backlog drain.

    HQ replays often disable live frames (``GRTECLYN_FRAMES=0``) to keep the
    watch consumer fast; PNG movies must still be rendered once from the
    remaining plotfiles.  Fast drains set ``GRTECLYN_CONSUMER_DRAIN=1``.
    """
    if explicit is not None:
        return explicit
    if _env_flag("GRTECLYN_CONSUMER_DRAIN"):
        return False
    return True


def consumer_delete_plotfiles_enabled(*, frames: bool, delete_requested: bool) -> bool:
    """Whether processed plotfiles may be deleted.

    Metrics-only live consumers (``frames=False``) must not delete AMR dumps
    unless ``GRTECLYN_DELETE_WITHOUT_FRAMES=1`` — otherwise post-run frame
    drains have nothing left to render.
    """
    if not delete_requested:
        return False
    if frames:
        return True
    return _env_flag("GRTECLYN_DELETE_WITHOUT_FRAMES")


def _resolve_n_points(default: int) -> int:
    """Angular resolution for spherical Psi4 extraction (HQ GW runs)."""
    if not _psi4_enabled():
        return default
    raw = os.environ.get("GRTECLYN_PSI4_N_POINTS", "").strip()
    if not raw:
        return 128
    try:
        return max(4, int(raw))
    except ValueError:
        return 128


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


def consumer_radii_from_env(
    default: Sequence[float] = (8.0, 12.0, 24.0),
) -> tuple[float, ...]:
    """Parse ``CONSUMER_RADII`` (space-separated) for HQ GW / shell extraction."""
    raw = os.environ.get("CONSUMER_RADII", "").split()
    if raw:
        return tuple(float(r) for r in raw)
    return tuple(float(r) for r in default)


def consumer_jobs_from_env(*, default: int = 1) -> int:
    """Parallel plotfile workers. Safe to raise on large-RAM nodes (``CONSUMER_JOBS=4``)."""
    raw = os.environ.get("CONSUMER_JOBS", "").strip()
    if not raw:
        jobs = default
    else:
        try:
            jobs = max(1, int(raw))
        except ValueError:
            jobs = default
    cap_raw = os.environ.get("GRTECLYN_CONSUMER_JOBS_MAX", "").strip()
    if cap_raw:
        try:
            jobs = min(jobs, max(1, int(cap_raw)))
        except ValueError:
            pass
    return jobs


def consumer_drain_minimal() -> bool:
    """Psi4 + delete only — skip FTL/geodesic/areal/shell/incremental score."""
    return _env_flag("GRTECLYN_CONSUMER_DRAIN_MINIMAL")


def build_consume_command(
    episode: Episode,
    *,
    profile: ConsumerProfile = "radial",
    radii: Sequence[float] = (8.0, 16.0),
    n_points: int = 64,
    delete: bool = True,
    keep_last: int = 1,
    watch: bool = True,
    jobs: int = 1,
    frames: bool = True,
    keep_existing_frames: bool = False,
    stable_seconds: float | None = None,
    ftl_timeseries: bool = False,
    central_timeseries: bool | None = None,
    confinement_timeseries: bool | None = None,
    ftl_L: float | None = None,
    incremental_score: bool = True,
    objective_mode: str = "weighted",
    target_stop_time: float | None = None,
    score_weights: Mapping[str, float] | None = None,
    evolving_geodesic: bool = False,
) -> list[str]:
    """Return argv for streaming plotfile extraction into episode/small_data."""
    center = _read_vector_param(episode.params_path, "center", (0.0, 0.0, 0.0))
    if len(center) < 3:
        center = (*center, *([0.0] * (3 - len(center))))
    center = center[:3]
    n_points = _resolve_n_points(n_points)
    l_full = _read_float_param(episode.params_path, "L_full", 40.0)
    # Frame/projection rendering is the dominant per-plotfile cost (12 yt renders
    # on a refined AMR grid, 10-70s each at max_level=3) and is *not* needed for
    # QD scoring -- only the FTL time-series + scalar metrics are.  Per-eval movies
    # are an inspection tool for promoted candidates, so the QD loop sets
    # GRTECLYN_FRAMES=0 to keep the consumer fast and stop the post-processing
    # backlog from starving the search.  Post-run drain sets GRTECLYN_CONSUMER_DRAIN=1.
    if not consumer_frames_enabled():
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
    ]
    if not consumer_drain_minimal():
        command.append("--areal-radius")
    command.extend(["-j", str(jobs), "--verbose"])

    if watch:
        command.append("--watch")
    effective_delete = consumer_delete_plotfiles_enabled(
        frames=frames, delete_requested=delete
    )
    if effective_delete:
        command.extend(["--delete", "--keep-last", str(keep_last)])
    if keep_existing_frames:
        command.append("--keep-existing-frames")
    if stable_seconds is not None:
        command.extend(["--stable-seconds", f"{float(stable_seconds):g}"])
    if ftl_timeseries:
        command.append("--ftl-timeseries")
        if ftl_L is not None:
            command.extend(["--ftl-l", f"{float(ftl_L):g}"])
        if evolving_geodesic:
            command.append("--evolving-geodesic")
    if _confinement_timeseries_enabled(confinement_timeseries, ftl_timeseries):
        well_width = _read_float_param(
            episode.params_path, "trajectory_well_width", 1.5
        )
        command.extend(
            ["--confinement-timeseries", "--confinement-well-width", f"{well_width:g}"]
        )
    central_enabled = _central_timeseries_enabled(central_timeseries)
    splash_incremental = (
        _incremental_score_enabled(incremental_score)
        and central_enabled
        and objective_mode == "critical_collapse"
    )
    if incremental_score and (ftl_timeseries or splash_incremental):
        stop_time = target_stop_time
        if stop_time is None:
            stop_time = _read_float_param(episode.params_path, "stop_time", 16.0)
        command.extend(
            [
                "--incremental-score",
                "--objective-mode",
                str(objective_mode),
                "--target-stop-time",
                f"{float(stop_time):g}",
            ]
        )
        if score_weights:
            for key, value in score_weights.items():
                command.extend(["--score-weight", f"{key}={value:g}"])
    if central_enabled:
        command.append("--central-timeseries")
        if _central_ball_enabled():
            command.append("--central-ball")
        if _central_radial_enabled():
            command.extend(["--central-radial-profile", "--central-radial-r-max", "6.0"])
        if _splash_early_term_enabled():
            command.extend(
                [
                    "--splash-early-term",
                    "--stop-sim-path",
                    str(episode.path / ".stop_sim"),
                ]
            )
        if not ftl_timeseries and objective_mode != "weighted":
            if "--objective-mode" not in command:
                command.extend(["--objective-mode", str(objective_mode)])

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
        if _psi4_enabled():
            command.append("--psi4")
        else:
            command.append("--no-psi4")
        if not consumer_drain_minimal():
            command.extend(
                [
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
            if _env_flag("GRTECLYN_FRAMES_AUTO_ZLIM"):
                command.extend(["--frames-auto-zlim", "--no-frames-global-zlim"])
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
