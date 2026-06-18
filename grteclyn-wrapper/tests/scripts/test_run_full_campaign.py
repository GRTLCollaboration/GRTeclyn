"""Tests for the end-to-end campaign orchestrator shell script."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path


WRAPPER_ROOT = Path(__file__).resolve().parents[2]
GRTECLYN_ROOT = WRAPPER_ROOT.parent
RUN_SCRIPT = WRAPPER_ROOT / "scripts/campaigns/run_full_campaign.sh"
PICK_SCRIPT = WRAPPER_ROOT / "scripts/campaigns/lib/pick_top_eval.py"


def _run_campaign(
    tmp_path: Path,
    *,
    campaign_name: str,
    extra_env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    campaign_root = tmp_path / "campaigns" / campaign_name
    env = {
        "GRTECLYN_ROOT": str(GRTECLYN_ROOT),
        "CAMPAIGN_ROOT": str(campaign_root),
        "CAMPAIGN_NAME": campaign_name,
        "DRY_RUN": "1",
        "QD_TARGET_EVALS": "1",
        "SKIP_QD_PREFLIGHT_TESTS": "1",
        "SKIP_CMA_PREFLIGHT_TESTS": "1",
        "GPU_IDS": "0",
        "BATCH_SIZE": "1",
    }
    if extra_env:
        env.update(extra_env)
    return subprocess.run(
        ["bash", str(RUN_SCRIPT)],
        cwd=WRAPPER_ROOT,
        env={**os.environ, **env},
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )


def test_run_full_campaign_script_syntax() -> None:
    subprocess.run(["bash", "-n", str(RUN_SCRIPT)], check=True)


def test_run_full_campaign_dry_run_qd_stage(tmp_path: Path) -> None:
    name = "pytest_qd_stage"
    result = _run_campaign(tmp_path, campaign_name=name, extra_env={"STAGE": "qd"})
    assert result.returncode == 0, result.stderr + result.stdout

    root = tmp_path / "campaigns" / name
    state = json.loads((root / "campaign_state.json").read_text(encoding="utf-8"))
    assert state["qd"]["status"] == "complete"
    assert (root / "campaign.json").is_file()
    assert "Stage 0: MAP-Elites" in result.stdout


def test_run_full_campaign_hq_dry_run(tmp_path: Path) -> None:
    name = "pytest_hq_stage"
    root = tmp_path / "campaigns" / name
    cmaes = root / "cmaes"
    cmaes.mkdir(parents=True)
    (cmaes / "trajectory.jsonl").write_text(
        '{"eval": 46, "status": "gpu_ok", "score": 179.8, "components": {}}\n',
        encoding="utf-8",
    )
    (cmaes / "eval_000046").mkdir()

    result = _run_campaign(
        tmp_path,
        campaign_name=name,
        extra_env={"STAGE": "hq"},
    )
    assert result.returncode == 0, result.stderr + result.stdout
    state = json.loads((root / "campaign_state.json").read_text(encoding="utf-8"))
    assert state["hq"]["status"] == "complete"
    assert state["hq"].get("dry_run") is True
    assert "Stage 2: HQ promotion" in result.stdout
    assert "Objective     : general_ftl" in result.stdout
    assert "Frames        : 1" in result.stdout


def test_run_full_campaign_resume_skips_qd(tmp_path: Path) -> None:
    name = "pytest_resume"
    root = tmp_path / "campaigns" / name
    root.mkdir(parents=True)
    (root / "logs").mkdir()
    (root / "campaign_state.json").write_text(
        json.dumps(
            {
                "qd": {"status": "complete", "top_eval": 63},
                "cmaes": {"status": "complete", "top_eval": 46},
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (root / "qd").mkdir()
    (root / "qd" / "trajectory.jsonl").write_text(
        '{"eval": 63, "status": "gpu_ok", "score": 165.0, "components": {}}\n',
        encoding="utf-8",
    )
    cmaes = root / "cmaes"
    cmaes.mkdir()
    (cmaes / "trajectory.jsonl").write_text(
        '{"eval": 46, "status": "gpu_ok", "score": 179.8, "components": {}}\n',
        encoding="utf-8",
    )
    (cmaes / "eval_000046").mkdir()

    result = _run_campaign(
        tmp_path,
        campaign_name=name,
        extra_env={"RESUME": "1", "STAGE": "all"},
    )
    assert result.returncode == 0, result.stderr + result.stdout
    assert "skip qd" in result.stderr
    assert "skip cmaes" in result.stderr
    assert "Stage 0: MAP-Elites" not in result.stdout
    assert "Stage 1: CMA-ES" not in result.stdout
    assert "Stage 2: HQ promotion" in result.stdout


def test_run_full_campaign_dry_run_all_stages(tmp_path: Path) -> None:
    name = "pytest_all_stages"
    result = _run_campaign(tmp_path, campaign_name=name)
    assert result.returncode == 0, result.stderr + result.stdout

    root = tmp_path / "campaigns" / name
    state = json.loads((root / "campaign_state.json").read_text(encoding="utf-8"))
    assert state["qd"]["status"] == "complete"
    assert state["cmaes"]["status"] == "complete"
    assert state["hq"]["status"] == "complete"
    assert "Stage 0: MAP-Elites" in result.stdout
    assert "Stage 1: CMA-ES" in result.stdout
    assert "Stage 2: HQ promotion" in result.stdout


def test_pick_top_eval_cli_pairs(tmp_path: Path) -> None:
    traj = tmp_path / "trajectory.jsonl"
    traj.write_text(
        '{"eval": 46, "status": "gpu_ok", "score": 179.8, "components": {}}\n',
        encoding="utf-8",
    )
    result = subprocess.run(
        ["python", str(PICK_SCRIPT), str(traj), "--format", "pairs", "--gpu", "0"],
        cwd=WRAPPER_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    assert result.stdout.strip() == "46 0"
