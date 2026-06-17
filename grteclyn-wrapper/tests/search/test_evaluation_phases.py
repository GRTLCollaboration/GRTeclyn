"""Tests for evaluate_overrides CPU/GPU phase split."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from grteclyn_wrapper.core.config import resolve_example
from grteclyn_wrapper.core.evaluation import (
    CpuGateResult,
    Evaluation,
    _run_cpu_grtresna_gates,
    _run_gpu_session,
    evaluate_overrides,
)
from grteclyn_wrapper.core.episode import create_episode
from grteclyn_wrapper.search.gpu_pool import GpuPool, GpuPoolConfig
from grteclyn_wrapper.search.grtresna_evaluation_gates import GRTresnaPreEvolutionGateConfig


@pytest.fixture
def episode(tmp_path: Path):
    return create_episode(tmp_path, name="phase_episode")


def test_cpu_gates_never_acquire_gpu(tmp_path: Path) -> None:
    pool = GpuPool(GpuPoolConfig(gpu_ids=[0], slots_per_gpu=1))
    with patch(
        "grteclyn_wrapper.core.evaluation._run_cpu_grtresna_gates",
        return_value=Evaluation(
            score=-80.0,
            components={"grtresna_rejection": -80.0},
            notes=["rejected"],
            episode_path=str(tmp_path),
            exit_code=None,
            preflight_rejected=False,
            reason="rejected",
            metrics={},
        ),
    ) as mock_cpu, patch.object(pool, "lease") as mock_lease:
        result = evaluate_overrides(
            {"grtresna_lumps": 1},
            out_dir=tmp_path,
            name="eval_000001",
            example=resolve_example("RadialRecipe"),
            template=resolve_example("RadialRecipe").template,
            executable=None,
            dry_run=True,
            grtresna=True,
            gpu_pool=pool,
        )
    assert result.reason == "rejected"
    mock_cpu.assert_called_once()
    mock_lease.assert_not_called()


def test_postload_reject_single_lease_released(tmp_path: Path, episode) -> None:
    gridinit = tmp_path / "initial_data.gridinit"
    gridinit.write_bytes(b"x")
    cpu_result = CpuGateResult(
        episode=episode,
        gte_overrides={"recipe_initial_data_file": str(gridinit)},
        gridinit_path=gridinit,
        gate_config=GRTresnaPreEvolutionGateConfig(postload_enabled=True),
    )
    pool = GpuPool(GpuPoolConfig(gpu_ids=[0], slots_per_gpu=1))
    acquires = 0
    original_lease = pool.lease

    def counting_lease():
        nonlocal acquires
        acquires += 1
        return original_lease()

    with patch.object(pool, "lease", side_effect=counting_lease), patch(
        "grteclyn_wrapper.core.evaluation.reject_postload_gate",
        return_value=MagicMock(
            reason="postload_fail",
            fitness=90.0,
            component_key="postload_rejection",
            metadata={"postload_rejected": True},
            notes=(),
            metrics=None,
        ),
    ), patch("grteclyn_wrapper.core.evaluation.run_episode") as mock_episode:
        with pool.lease() as gpu:
            result = _run_gpu_session(
                cpu_result=cpu_result,
                example=resolve_example("RadialRecipe"),
                template=resolve_example("RadialRecipe").template,
                executable=None,
                cuda_devices=gpu,
                check_params=False,
                dry_run=True,
                target_stop_time=1.0,
                score_weights=None,
                objective_mode="weighted",
                ftl_L=None,
                consume_plotfiles=False,
                consumer_radii=(4.0,),
                consumer_keep_last=1,
                consumer_ftl_timeseries=False,
                consumer_evolving_geodesic=None,
                grtresna_postload_gate=True,
                postload_gate_config=None,
            )
    assert result.reason == "postload_fail"
    mock_episode.assert_not_called()
    assert acquires == 1


def test_postload_pass_lease_spans_evolution(tmp_path: Path, episode) -> None:
    gridinit = tmp_path / "initial_data.gridinit"
    gridinit.write_bytes(b"x")
    cpu_result = CpuGateResult(
        episode=episode,
        gte_overrides={"recipe_initial_data_file": str(gridinit)},
        gridinit_path=gridinit,
        gate_config=GRTresnaPreEvolutionGateConfig(postload_enabled=True),
    )
    pool = GpuPool(GpuPoolConfig(gpu_ids=[0], slots_per_gpu=1))
    lease_count = 0
    original_lease = pool.lease

    def counting_lease():
        nonlocal lease_count
        lease_count += 1
        return original_lease()

    with patch.object(pool, "lease", side_effect=counting_lease), patch(
        "grteclyn_wrapper.core.evaluation.reject_postload_gate",
        return_value=None,
    ), patch(
        "grteclyn_wrapper.core.evaluation.run_episode",
        return_value=MagicMock(returncode=0),
    ) as mock_episode, patch(
        "grteclyn_wrapper.core.evaluation.read_episode_metrics",
        return_value=MagicMock(__dataclass_fields__={}),
    ), patch(
        "grteclyn_wrapper.core.evaluation.score_episode",
        return_value=MagicMock(total=2.0, components={}, notes=[]),
    ), patch(
        "grteclyn_wrapper.core.evaluation.dataclass_to_dict",
        side_effect=lambda obj: {"total": 2.0, "components": {}, "notes": []}
        if hasattr(obj, "total")
        else {},
    ), patch(
        "grteclyn_wrapper.core.evaluation.write_params",
    ):
        with pool.lease() as gpu:
            _run_gpu_session(
                cpu_result=cpu_result,
                example=resolve_example("RadialRecipe"),
                template=resolve_example("RadialRecipe").template,
                executable=MagicMock(),
                cuda_devices=gpu,
                check_params=False,
                dry_run=False,
                target_stop_time=1.0,
                score_weights=None,
                objective_mode="weighted",
                ftl_L=None,
                consume_plotfiles=False,
                consumer_radii=(4.0,),
                consumer_keep_last=1,
                consumer_ftl_timeseries=False,
                consumer_evolving_geodesic=None,
                grtresna_postload_gate=True,
                postload_gate_config=None,
            )
    assert lease_count == 1
    mock_episode.assert_called_once()


def test_postload_disabled_skips_gpu_gate(tmp_path: Path, episode) -> None:
    gridinit = tmp_path / "initial_data.gridinit"
    gridinit.write_bytes(b"x")
    cpu_result = CpuGateResult(
        episode=episode,
        gte_overrides={"recipe_initial_data_file": str(gridinit)},
        gridinit_path=gridinit,
        gate_config=GRTresnaPreEvolutionGateConfig(),
    )
    with patch(
        "grteclyn_wrapper.core.evaluation.reject_postload_gate",
    ) as mock_postload, patch(
        "grteclyn_wrapper.core.evaluation.read_episode_metrics",
        return_value=MagicMock(__dataclass_fields__={}),
    ), patch(
        "grteclyn_wrapper.core.evaluation.score_episode",
        return_value=MagicMock(total=1.0, components={}, notes=[]),
    ), patch(
        "grteclyn_wrapper.core.evaluation.dataclass_to_dict",
        side_effect=lambda obj: {"total": 1.0, "components": {}, "notes": []}
        if hasattr(obj, "total")
        else {},
    ), patch(
        "grteclyn_wrapper.core.evaluation.write_params",
    ):
        _run_gpu_session(
            cpu_result=cpu_result,
            example=resolve_example("RadialRecipe"),
            template=resolve_example("RadialRecipe").template,
            executable=None,
            cuda_devices="0",
            check_params=False,
            dry_run=True,
            target_stop_time=1.0,
            score_weights=None,
            objective_mode="weighted",
            ftl_L=None,
            consume_plotfiles=False,
            consumer_radii=(4.0,),
            consumer_keep_last=1,
            consumer_ftl_timeseries=False,
            consumer_evolving_geodesic=None,
            grtresna_postload_gate=False,
            postload_gate_config=None,
        )
    mock_postload.assert_not_called()
