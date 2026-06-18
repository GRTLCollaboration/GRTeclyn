"""Tests for pre-GPU rejection learning (QD near-miss pool, shadow archive, CMA warm-start)."""

from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from grteclyn_wrapper.search.grtresna_convergence_gate import GRTresnaConvergenceConfig
from grteclyn_wrapper.search.grtresna_evaluation_gates import reject_grtresna_convergence
from grteclyn_wrapper.search.optimize import SearchDimension
from grteclyn_wrapper.search.optimize.candidates import _load_warm_start_vectors
from grteclyn_wrapper.search.pre_gpu import (
    NearMissPool,
    is_graded_rejection,
    parse_convergence_from_reason,
    pre_gpu_descriptor_details,
)
from grteclyn_wrapper.search.qd_search.archive import QDArchive
from grteclyn_wrapper.search.qd_search.pre_gpu_archive import rebuild_pre_gpu_archive_from_trajectory
from grteclyn_wrapper.search.qd_search.sampling import QdSamplingConfig, sample_next_candidate


def _dim(key: str = "grtresna_shell_amp") -> list[SearchDimension]:
    return [SearchDimension(param_key=key, lower=0.0, upper=1.0, initial=0.5)]


def test_reject_grtresna_convergence_includes_metrics() -> None:
    rejection = reject_grtresna_convergence(
        {"ham_pct": 8.0, "mom_pct": 3.0},
        config=GRTresnaConvergenceConfig(),
    )
    assert rejection is not None
    assert rejection.metrics is not None
    assert rejection.metrics["grtresna_convergence"]["ham_pct"] == 8.0


def test_parse_convergence_from_reason_legacy_v21() -> None:
    parsed = parse_convergence_from_reason(
        "GRTresna convergence too poor: Ham=100% Mom=0% thresholds=(5%, 5%)"
    )
    assert parsed == {"ham_pct": 100.0, "mom_pct": 0.0}


def test_pre_gpu_descriptors_ham_mom_axes_differ() -> None:
    mild = pre_gpu_descriptor_details(
        {
            "status": "grtresna_rejected",
            "reason": "Ham=6% Mom=4%",
            "score": -110.0,
            "overrides": {},
        },
        convergence_config=GRTresnaConvergenceConfig(max_ham_pct=5.0, max_mom_pct=5.0),
    )
    severe = pre_gpu_descriptor_details(
        {
            "status": "grtresna_rejected",
            "reason": "Ham=100% Mom=80%",
            "score": -350.0,
            "overrides": {},
        },
        convergence_config=GRTresnaConvergenceConfig(max_ham_pct=5.0, max_mom_pct=5.0),
    )
    assert mild["x"] > severe["x"]
    assert mild["y"] > severe["y"]


def test_near_miss_pool_keeps_best_and_excludes_failed() -> None:
    pool = NearMissPool(max_size=2)
    pool.consider({
        "status": "grtresna_rejected",
        "score": -350.0,
        "overrides": {"grtresna_shell_amp": 0.1},
    })
    pool.consider({
        "status": "grtresna_rejected",
        "score": -195.0,
        "overrides": {"grtresna_shell_amp": 0.2},
    })
    pool.consider({
        "status": "grtresna_failed",
        "score": -350.0,
        "overrides": {"grtresna_shell_amp": 0.9},
    })
    assert len(pool) == 2
    params = pool.param_sets()
    assert all(p["grtresna_shell_amp"] != 0.9 for p in params)
    assert any(abs(p["grtresna_shell_amp"] - 0.2) < 1e-9 for p in params)


def test_shadow_archive_distinct_cells_for_ham_mom() -> None:
    records = [
        {
            "status": "grtresna_rejected",
            "score": -120.0,
            "reason": "Ham=6% Mom=4%",
            "overrides": {"grtresna_shell_amp": 0.1},
        },
        {
            "status": "grtresna_rejected",
            "score": -300.0,
            "reason": "Ham=100% Mom=80%",
            "overrides": {"grtresna_shell_amp": 0.2},
        },
    ]
    archive = rebuild_pre_gpu_archive_from_trajectory(records, bins=8)
    cells = {e.cell for e in archive.cells.values()}
    assert len(cells) >= 2


def test_sample_next_candidate_near_miss_after_rejections() -> None:
    rng = np.random.default_rng(0)
    archive = QDArchive(bins=8)
    shadow = rebuild_pre_gpu_archive_from_trajectory(
        [
            {
                "status": "grtresna_rejected",
                "score": -195.0,
                "reason": "Ham=100% Mom=0%",
                "overrides": {"grtresna_shell_amp": 0.42},
            }
        ],
        bins=8,
    )
    pool = NearMissPool.rebuild_from_trajectory(
        [
            {
                "status": "grtresna_rejected",
                "score": -195.0,
                "overrides": {"grtresna_shell_amp": 0.42},
            }
        ],
        max_size=8,
    )
    cfg = QdSamplingConfig(
        elite_fraction=0.0,
        shadow_fraction=0.0,
        near_miss_fraction=1.0,
        feasible_fraction=0.0,
        random_fraction=0.0,
        mutation_sigma=0.01,
    )
    child = sample_next_candidate(
        dims=_dim(),
        rng=rng,
        archive=archive,
        pre_gpu_archive=shadow,
        near_miss_pool=pool,
        config=cfg,
    )
    assert abs(child[0] - 0.42) < 0.05


def test_load_warm_start_includes_near_miss_vectors() -> None:
    dims = _dim()
    with TemporaryDirectory() as tmp:
        traj = Path(tmp) / "trajectory.jsonl"
        traj.write_text(
            json.dumps({
                "status": "gpu_ok",
                "score": 10.0,
                "overrides": {"grtresna_shell_amp": 0.11},
            }) + "\n"
            + json.dumps({
                "status": "grtresna_rejected",
                "score": -195.0,
                "overrides": {"grtresna_shell_amp": 0.42},
            }) + "\n",
            encoding="utf-8",
        )
        vectors = _load_warm_start_vectors(
            [traj], dims, top_k=1, include_near_miss=True, near_miss_k=4,
        )
    assert len(vectors) == 2
    amps = sorted(vec[0] for vec in vectors)
    assert amps[0] == 0.11
    assert amps[1] == 0.42


def test_is_graded_rejection_classification() -> None:
    assert is_graded_rejection({"status": "grtresna_rejected"})
    assert not is_graded_rejection({"status": "grtresna_failed"})


def test_hq_replay_does_not_import_pre_gpu_learning() -> None:
    replay_path = (
        Path(__file__).resolve().parents[2]
        / "scripts/campaigns/hq/replay_eval.py"
    )
    source = replay_path.read_text(encoding="utf-8")
    assert "pre_gpu" not in source
    assert "near_miss_pool" not in source
