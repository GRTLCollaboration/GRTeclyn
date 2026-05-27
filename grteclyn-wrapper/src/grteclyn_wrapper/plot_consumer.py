"""Build consume_plotfiles command lines for GRTeclyn wrapper episodes."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path
from typing import Literal, Sequence

from .config import REPO_ROOT
from .episode import Episode

ConsumerProfile = Literal["wormhole", "radial"]


def resolve_consume_python() -> list[str]:
    """Return argv prefix for Python with visualization deps (yt)."""
    if shutil.which("uv") is not None and (REPO_ROOT / "pyproject.toml").exists():
        return ["uv", "run", "python"]
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
) -> list[str]:
    """Return argv for streaming plotfile extraction into episode/small_data."""
    command = [
        *resolve_consume_python(),
        "-m",
        "src.visualisation.process_wave.consume_plotfiles",
        "--data",
        str(episode.path),
        "--out",
        str(episode.small_data_dir),
        "--radii",
        *[str(r) for r in radii],
        "--n-points",
        str(n_points),
        "--center",
        "0.0",
        "0.0",
        "0.0",
        "--areal-radius",
        "-j",
        str(jobs),
        "--verbose",
    ]

    if watch:
        command.append("--watch")
    if delete:
        command.extend(["--delete", "--keep-last", str(keep_last)])

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
            command.extend(
                [
                    "--frames-fields",
                    "chi",
                    "lapse",
                    "K",
                    "--frames-axis",
                    "z",
                    "--frames-out",
                    str(episode.frames_dir),
                ]
            )

    return command


def default_radii_for_example(example_name: str) -> tuple[float, ...]:
    if example_name == "RadialRecipe":
        return (4.0, 8.0)
    return (12.0, 16.0, 20.0, 24.0)
