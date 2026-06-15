from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping


def bounded_reward(value: float, scale: float) -> float:
    if not math.isfinite(value) or value < 0:
        return 0.0
    return 1.0 / (1.0 + value / scale)


def graded(x: float, scale: float, ceiling: float = 1.0) -> float:
    if not math.isfinite(x) or x <= 0.0:
        return 0.0
    return min(math.log1p(x / scale), ceiling)


def domain_half_width_from_overrides(
    overrides: Mapping[str, object] | None,
) -> float | None:
    """Half-width of the cubic domain (``L_full / 2``) for the horizon-finder
    sanity guard.  The domain is ``[0, L_full]`` with the geometry centered at
    ``L_full / 2``, so a genuine interior trapped surface must satisfy
    ``r << L_full / 2``."""
    if not overrides:
        return None
    try:
        l_full = float(overrides.get("L_full"))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return 0.5 * l_full if l_full > 0.0 else None


def domain_half_width_from_params(params_path: Path) -> float | None:
    """Read ``L_full`` from a written ``params.txt`` when overrides omit it."""
    if not params_path.is_file():
        return None
    for line in params_path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped.startswith("L_full"):
            continue
        if "=" not in stripped:
            continue
        try:
            l_full = float(stripped.split("=", 1)[1].strip())
        except ValueError:
            return None
        return 0.5 * l_full if l_full > 0.0 else None
    return None


def domain_half_width_for_episode(
    episode_dir: Path,
    overrides: Mapping[str, object] | None = None,
) -> float | None:
    """Resolve domain half-width from overrides or the episode ``params.txt``."""
    half_width = domain_half_width_from_overrides(overrides)
    if half_width is not None:
        return half_width
    return domain_half_width_from_params(episode_dir / "params.txt")
