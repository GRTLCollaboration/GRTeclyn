from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from .config import REPO_ROOT


def _git_commit() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return None
    return result.stdout.strip() or None


def _next_episode_dir(runs_dir: Path) -> Path:
    runs_dir.mkdir(parents=True, exist_ok=True)
    used = []
    for child in runs_dir.iterdir():
        if child.is_dir() and child.name.startswith("episode_"):
            try:
                used.append(int(child.name.removeprefix("episode_")))
            except ValueError:
                continue
    return runs_dir / f"episode_{(max(used) + 1) if used else 1:06d}"


@dataclass(frozen=True)
class Episode:
    path: Path

    @property
    def params_path(self) -> Path:
        return self.path / "params.txt"

    @property
    def metadata_path(self) -> Path:
        return self.path / "metadata.json"

    @property
    def log_path(self) -> Path:
        return self.path / "run.log"

    @property
    def score_path(self) -> Path:
        return self.path / "score.json"

    @property
    def data_dir(self) -> Path:
        return self.path / "data"

    @property
    def small_data_dir(self) -> Path:
        return self.path / "small_data"

    @property
    def frames_dir(self) -> Path:
        return self.path / "frames"


def create_episode(
    runs_dir: Path,
    *,
    name: str | None = None,
    metadata: Mapping[str, Any] | None = None,
) -> Episode:
    episode_dir = (runs_dir / name) if name else _next_episode_dir(runs_dir)
    episode_dir = episode_dir.expanduser().resolve()
    episode_dir.mkdir(parents=True, exist_ok=False)
    (episode_dir / "data").mkdir()
    (episode_dir / "small_data").mkdir()

    payload: dict[str, Any] = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(),
        "repo_root": str(REPO_ROOT),
    }
    if metadata:
        payload.update(metadata)
    write_json(episode_dir / "metadata.json", payload)
    return Episode(episode_dir)


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def update_metadata(episode: Episode, values: Mapping[str, Any]) -> None:
    if episode.metadata_path.exists():
        payload = json.loads(episode.metadata_path.read_text(encoding="utf-8"))
    else:
        payload = {}
    payload.update(values)
    write_json(episode.metadata_path, payload)
