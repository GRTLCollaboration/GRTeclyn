#!/usr/bin/env python3
"""Enforce the pump convention on a promote-campaign manifest.

Convention (research/neuralspacetime/GPU_RUN_PLAN.md section 12.1): the PD
pump runs for the ENTIRE simulation, never stopped mid-run. In config terms
that means NO ``rl_pump_stop_time`` anywhere -- an absent key IS the
convention: the wiring writes the key only for values >= 0, and the
evolution's built-in default is -1 = never stop.

A manifest that deliberately stops the pump (a pump-off control such as the
pump-free twin) must say so out loud with a top-level

    "pump_off_control": true

and only then are non-negative ``rl_pump_stop_time`` values accepted.  The
flag is intentionally NOT inherited silently: a clone that copies stale
pump values without copying the declaration is refused at launch.

Pump-on manifests must additionally have the emission floor pinned in the
environment (``GEODESIC_EMIT_MIN_TIME``): with no non-negative pump value the
scorer's fallback floor vanishes, every launch from t=0 counts, and f_geo
silently changes meaning (metrics/score/ftl.py).

Exit codes: 0 ok, 3 refused, 2 usage.  Wired into run_matrix.sh for every
--list and launch; standalone:

    python3 validate_pump_convention.py <manifest.json>
"""

from __future__ import annotations

import json
import math
import os
import sys
from typing import Any, Iterator

PUMP_KEY = "rl_pump_stop_time"
FLOOR_ENV = "GEODESIC_EMIT_MIN_TIME"
PUMP_ENV = "RL_PUMP_STOP_TIME"


def find_pump_settings(node: Any, path: str = "manifest") -> Iterator[tuple[str, float]]:
    """Yield (json_path, value) for every rl_pump_stop_time in the tree.

    Catches both the dict form ({"rl_pump_stop_time": 4.0}, e.g. in
    physics_frozen) and the override-string form ("rl_pump_stop_time=4",
    e.g. in extra_overrides).  Unparseable values yield NaN, which the
    caller treats as a violation rather than a pass.
    """
    if isinstance(node, dict):
        for key, value in node.items():
            child = f"{path}.{key}"
            if key == PUMP_KEY:
                try:
                    yield child, float(value)
                except (TypeError, ValueError):
                    yield child, math.nan
            else:
                yield from find_pump_settings(value, child)
    elif isinstance(node, list):
        for index, value in enumerate(node):
            yield from find_pump_settings(value, f"{path}[{index}]")
    elif isinstance(node, str) and f"{PUMP_KEY}=" in node:
        raw = node.split(f"{PUMP_KEY}=", 1)[1].split()[0]
        try:
            yield f"{path} ({node!r})", float(raw)
        except ValueError:
            yield f"{path} ({node!r})", math.nan


def validate(manifest_path: str, env: dict[str, str]) -> list[str]:
    """Return a list of violations (empty = compliant)."""
    with open(manifest_path, encoding="utf-8") as fh:
        manifest = json.load(fh)

    # NaN must not pass: "not (v < 0)" is True for NaN.
    stops = [(p, v) for p, v in find_pump_settings(manifest) if not v < 0.0]
    opt_in = manifest.get("pump_off_control") is True
    errors: list[str] = []

    if stops and not opt_in:
        listing = "\n".join(f"      {p} -> {v}" for p, v in stops)
        errors.append(
            "manifest stops the pump mid-run:\n"
            f"{listing}\n"
            "    The pump runs for the ENTIRE simulation (GPU_RUN_PLAN.md 12.1).\n"
            f"    Remove {PUMP_KEY} everywhere (absent key = pump always on), or\n"
            "    declare a deliberate pump-off control with top-level\n"
            '    "pump_off_control": true.'
        )

    env_raw = env.get(PUMP_ENV, "").strip()
    if env_raw and not opt_in:
        try:
            env_stop = float(env_raw)
        except ValueError:
            errors.append(f"env {PUMP_ENV}={env_raw!r} is not a number")
        else:
            if env_stop >= 0.0:
                errors.append(
                    f"env {PUMP_ENV}={env_raw} stops the pump mid-run; use -1\n"
                    '    (or declare "pump_off_control": true in the manifest).'
                )

    pump_on = not stops
    if pump_on and not env.get(FLOOR_ENV, "").strip():
        errors.append(
            f"pump-on manifest without {FLOOR_ENV} in the environment:\n"
            "    with no non-negative pump value the scorer's fallback emission\n"
            "    floor vanishes and every launch from t=0 counts -- f_geo silently\n"
            f"    changes meaning (metrics/score/ftl.py). Export {FLOOR_ENV}=4\n"
            "    in campaign.env.sh."
        )

    return errors


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(f"usage: {argv[0]} <manifest.json>", file=sys.stderr)
        return 2
    errors = validate(argv[1], dict(os.environ))
    if errors:
        print(f"[pump-convention] REFUSED {argv[1]}:", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
