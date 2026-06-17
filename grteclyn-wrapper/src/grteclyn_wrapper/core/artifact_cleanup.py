"""Tiered episode artifact cleanup to reclaim disk after scoring."""

from __future__ import annotations

import json
import logging
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

logger = logging.getLogger(__name__)

Tier = Literal["plotfiles_only", "full_non_hq"]

_PLOTFILE_PREFIXES = ("RadialRecipePlt", "RadialRecipeChk")
_GRIDINIT_NAMES = ("initial_data.gridinit", "initial_data.matter.json")
_GRTRESNA_DIR = "grtresna"
_PROTECTED_FILES = {"score.json", "metadata.json", "params.txt"}


@dataclass
class CleanupReport:
    episode_dir: Path
    tier: Tier
    removed_paths: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def _read_score_general_ftl_solved(score_path: Path) -> object | None:
    if not score_path.is_file():
        return None
    try:
        payload = json.loads(score_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    score = payload.get("score")
    if isinstance(score, dict):
        return score.get("general_ftl_solved")
    return None


def cleanup_episode_artifacts(
    episode_dir: Path,
    *,
    tier: Tier,
    keep_for_hq: bool = False,
) -> CleanupReport:
    """Remove heavy artifacts from a scored episode directory."""
    report = CleanupReport(episode_dir=episode_dir, tier=tier)
    if not episode_dir.is_dir():
        return report

    score_path = episode_dir / "score.json"
    if tier == "plotfiles_only" and not score_path.is_file():
        report.warnings.append("skipped plotfile cleanup: score.json missing")
        return report

    for child in episode_dir.iterdir():
        if child.is_dir() and any(child.name.startswith(prefix) for prefix in _PLOTFILE_PREFIXES):
            shutil.rmtree(child, ignore_errors=True)
            report.removed_paths.append(str(child))

    if tier == "plotfiles_only":
        return report

    if keep_for_hq:
        return report

    solved = _read_score_general_ftl_solved(score_path)
    for name in _GRIDINIT_NAMES:
        path = episode_dir / name
        if path.is_file():
            if solved is not None:
                report.warnings.append(
                    f"removing gridinit for eval with general_ftl_solved in score.json: {path}"
                )
                logger.warning(report.warnings[-1])
            path.unlink(missing_ok=True)
            report.removed_paths.append(str(path))

    grtresna_dir = episode_dir / _GRTRESNA_DIR
    if grtresna_dir.is_dir():
        shutil.rmtree(grtresna_dir, ignore_errors=True)
        report.removed_paths.append(str(grtresna_dir))

    return report
