"""Tests for the promote-manifest pump-convention gate.

Convention under test (GPU_RUN_PLAN.md section 12.1): the pump runs for the
entire simulation; an absent rl_pump_stop_time IS the convention; a manifest
that stops the pump must declare "pump_off_control": true; a pump-on
manifest requires GEODESIC_EMIT_MIN_TIME in the environment.
"""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

WRAPPER_ROOT = Path(__file__).resolve().parents[2]
VALIDATOR = (
    WRAPPER_ROOT / "scripts/campaigns/promote/lib/validate_pump_convention.py"
)
TEMPLATE_DIR = WRAPPER_ROOT / "scripts/campaigns/promote/bicomplex_cmaes_v1"

_spec = importlib.util.spec_from_file_location("validate_pump_convention", VALIDATOR)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

FLOOR_ENV = {"GEODESIC_EMIT_MIN_TIME": "4"}


def _write(tmp_path: Path, manifest: dict) -> str:
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    return str(path)


def test_baked_stop_in_physics_frozen_refused(tmp_path):
    path = _write(tmp_path, {"physics_frozen": {"rl_pump_stop_time": 4.0}})
    errors = _mod.validate(path, dict(FLOOR_ENV))
    assert any("stops the pump" in e for e in errors)


def test_baked_stop_in_extra_overrides_refused(tmp_path):
    path = _write(tmp_path, {"extra_overrides": ["rl_pump_stop_time=4"]})
    errors = _mod.validate(path, dict(FLOOR_ENV))
    assert any("stops the pump" in e for e in errors)


def test_unparseable_stop_refused_not_passed(tmp_path):
    path = _write(tmp_path, {"extra_overrides": ["rl_pump_stop_time=abc"]})
    errors = _mod.validate(path, dict(FLOOR_ENV))
    assert any("stops the pump" in e for e in errors)


def test_declared_pump_off_control_accepted(tmp_path):
    path = _write(
        tmp_path,
        {
            "pump_off_control": True,
            "physics_frozen": {"rl_pump_stop_time": 0.0},
            "extra_overrides": ["rl_pump_stop_time=0"],
        },
    )
    assert _mod.validate(path, {}) == []


def test_negative_stop_is_the_convention(tmp_path):
    path = _write(tmp_path, {"extra_overrides": ["rl_pump_stop_time=-1"]})
    assert _mod.validate(path, dict(FLOOR_ENV)) == []


def test_clean_manifest_with_floor_accepted(tmp_path):
    path = _write(tmp_path, {"extra_overrides": ["grtresna_scalar_mass=1.0"]})
    assert _mod.validate(path, dict(FLOOR_ENV)) == []


def test_clean_manifest_without_floor_refused(tmp_path):
    path = _write(tmp_path, {})
    errors = _mod.validate(path, {})
    assert any("GEODESIC_EMIT_MIN_TIME" in e for e in errors)


def test_env_stop_refused_without_declaration(tmp_path):
    path = _write(tmp_path, {})
    env = {"RL_PUMP_STOP_TIME": "4", **FLOOR_ENV}
    errors = _mod.validate(path, env)
    assert any("RL_PUMP_STOP_TIME=4" in e for e in errors)


def test_env_never_stop_accepted(tmp_path):
    path = _write(tmp_path, {})
    env = {"RL_PUMP_STOP_TIME": "-1", **FLOOR_ENV}
    assert _mod.validate(path, env) == []


def test_template_manifest_refused_end_to_end():
    """The retired bicomplex manifest bakes rl_pump_stop_time=4 twice and
    must be refused at launch instead of silently stopping the pump."""
    proc = subprocess.run(
        [sys.executable, str(VALIDATOR), str(TEMPLATE_DIR / "manifest.json")],
        env=FLOOR_ENV,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 3
    assert "REFUSED" in proc.stderr


def test_pumpfree_twin_accepted_end_to_end():
    proc = subprocess.run(
        [sys.executable, str(VALIDATOR), str(TEMPLATE_DIR / "manifest_pumpfree.json")],
        env={},
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
