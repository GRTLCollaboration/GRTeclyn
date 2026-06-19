"""Tests for splash / critical_collapse campaign shell wiring."""

from __future__ import annotations

import subprocess
from pathlib import Path


WRAPPER_ROOT = Path(__file__).resolve().parents[2]
SEARCH_COMMON = WRAPPER_ROOT / "scripts/campaigns/lib/search_common.sh"
SPLASH_RUN = WRAPPER_ROOT / "scripts/campaigns/splash/run.sh"


def _source_env(script: str, extra: str = "") -> dict[str, str]:
    cmd = f"""
set -euo pipefail
{extra}
source "{script}"
printf '%s\\n' "OBJECTIVE_MODE=${{OBJECTIVE_MODE:-}}" \
  "GRTECLYN_EVOLVING_GEODESIC=${{GRTECLYN_EVOLVING_GEODESIC:-}}" \
  "GRTECLYN_CENTRAL_TIMESERIES=${{GRTECLYN_CENTRAL_TIMESERIES:-}}" \
  "GRTRESNA_MATTER_SECTOR=${{GRTRESNA_MATTER_SECTOR:-}}" \
  "GRTRESNA_MATTER_COUPLING=${{GRTRESNA_MATTER_COUPLING:-}}" \
  "GRTRESNA_ANSATZ=${{GRTRESNA_ANSATZ:-}}" \
  "DESCRIPTOR_MODE=${{DESCRIPTOR_MODE:-}}" \
  "SPLASH_MODE=${{SPLASH_MODE:-}}"
"""
    result = subprocess.run(
        ["bash", "-c", cmd],
        cwd=WRAPPER_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    env: dict[str, str] = {}
    for line in result.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            env[key] = value
    return env


def test_search_common_critical_collapse_flags() -> None:
    env = _source_env(str(SEARCH_COMMON), "export OBJECTIVE_MODE=critical_collapse")
    assert env["GRTECLYN_EVOLVING_GEODESIC"] == "0"
    assert env["GRTECLYN_CENTRAL_TIMESERIES"] == "1"


def test_search_common_ftl_first_unaffected() -> None:
    env = _source_env(str(SEARCH_COMMON), "export OBJECTIVE_MODE=ftl_first")
    assert env["GRTECLYN_EVOLVING_GEODESIC"] == "1"
    assert env.get("GRTECLYN_CENTRAL_TIMESERIES", "") != "1"


def test_search_common_exotic_coupling_allowed_outside_splash() -> None:
    env = _source_env(
        str(SEARCH_COMMON),
        "export OBJECTIVE_MODE=ftl_first\nexport GRTRESNA_MATTER_COUPLING=exotic",
    )
    assert env["GRTRESNA_MATTER_COUPLING"] == "exotic"


def test_splash_run_script_defaults() -> None:
    subprocess.run(["bash", "-n", str(SPLASH_RUN)], check=True)
    text = SPLASH_RUN.read_text(encoding="utf-8")
    assert "GRTRESNA_MATTER_SECTOR=scalar" in text
    assert 'OBJECTIVE_MODE="${OBJECTIVE_MODE:-critical_collapse}"' in text
    assert "GRTRESNA_ANSATZ=shell" in text
    assert 'DESCRIPTOR_MODE="${DESCRIPTOR_MODE:-wave_focusing}"' in text
    assert "export GRTECLYN_FRAMES=1" in text
    assert "GRTRESNA_MATTER_COUPLING=canonical" in text
    assert "GRTECLYN_FRAMES_AUTO_ZLIM=1" in text


def test_parser_accepts_splash_cli_choices() -> None:
    from grteclyn_wrapper.cli.parser import build_parser

    parser = build_parser()
    args = parser.parse_args(
        [
            "qd",
            "--grtresna",
            "--objective-mode",
            "critical_collapse",
            "--descriptor-mode",
            "wave_focusing",
            "--grtresna-ansatz",
            "splash",
        ]
    )
    assert args.objective_mode == "critical_collapse"
    assert args.descriptor_mode == "wave_focusing"
    assert args.grtresna_ansatz == "splash"
