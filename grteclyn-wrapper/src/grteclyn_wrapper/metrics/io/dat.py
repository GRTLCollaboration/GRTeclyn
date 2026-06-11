"""Parsing helpers for C++ diagnostic ``.dat`` files."""

from __future__ import annotations

from pathlib import Path


def numeric_rows(path: Path, min_columns: int) -> list[list[float]]:
    rows: list[list[float]] = []
    if not path.exists():
        return rows

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            try:
                values = [float(item) for item in stripped.split()]
            except ValueError:
                continue
            if len(values) >= min_columns:
                rows.append(values)
    rows.sort(key=lambda row: row[0])
    return rows


def resolve_episode_dat_path(episode_dir: Path, *relative_parts: str) -> Path:
    """Resolve a diagnostic path under ``data/`` with episode-root fallback."""
    data_path = episode_dir.joinpath("data", *relative_parts)
    if data_path.exists():
        return data_path
    return episode_dir.joinpath(*relative_parts)
