from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_elite_overrides(eval_dir: Path) -> dict[str, Any]:
    metadata_path = eval_dir / "metadata.json"
    if not metadata_path.exists():
        raise FileNotFoundError(f"Missing metadata.json in {eval_dir}")
    payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    overrides = payload.get("overrides")
    if not isinstance(overrides, dict):
        raise ValueError(f"metadata.json in {eval_dir} has no overrides dict")
    return dict(overrides)


def pump_geometry_from_overrides(overrides: dict[str, Any]) -> tuple[float, float]:
    radius = float(overrides.get("recipe_basis_radius_max", 5.0))
    width = float(overrides.get("recipe_basis_width", 1.5))
    return radius, width
