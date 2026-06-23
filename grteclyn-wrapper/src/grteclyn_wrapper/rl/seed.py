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


def _read_matter_metadata(eval_dir: Path) -> dict[str, Any]:
    path = eval_dir / "initial_data.matter.json"
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    return payload if isinstance(payload, dict) else {}


def lump_seeds_from_eval(
    eval_dir: Path, *, max_lumps: int = 5
) -> tuple[int, list[tuple[float, float, float]]]:
    """Return ``(num_lumps, centres)`` to seed the RL lump tracker.

    Centres are throat-relative offsets, i.e. positions in the RadialRecipe
    centred coordinate frame (the same frame the pump and tracker use).  Falls
    back to a single spotlight at the throat centre when no per-lump centres
    are recorded (e.g. a single coherent boson star).
    """
    matter = _read_matter_metadata(eval_dir)
    centers: list[tuple[float, float, float]] = []
    for c in matter.get("lump_centers") or ():
        if c is None or len(c) < 3:
            continue
        centers.append((float(c[0]), float(c[1]), float(c[2])))
        if len(centers) >= max_lumps:
            break
    if not centers:
        return 1, [(0.0, 0.0, 0.0)]
    return len(centers), centers


def rl_pump_params(
    eval_dir: Path,
    *,
    pump_width: float | None = None,
    max_amplitude: float = 0.05,
    max_lumps: int = 5,
) -> dict[str, Any]:
    """GRTeclyn ``params.txt`` entries configuring the RL multi-site pump.

    Lump seed positions / count come from the elite matter metadata; the
    spotlight envelope width defaults to the boson profile width when present.
    """
    num_lumps, centers = lump_seeds_from_eval(eval_dir, max_lumps=max_lumps)
    if pump_width is None:
        matter = _read_matter_metadata(eval_dir)
        pump_width = float(matter.get("bs_profile_width") or 0.0) or 1.5
    return {
        "rl_num_lumps": num_lumps,
        "rl_lump_seed_x": " ".join(repr(c[0]) for c in centers),
        "rl_lump_seed_y": " ".join(repr(c[1]) for c in centers),
        "rl_lump_seed_z": " ".join(repr(c[2]) for c in centers),
        "rl_pump_width": pump_width,
        "rl_pump_max_amplitude": max_amplitude,
    }
