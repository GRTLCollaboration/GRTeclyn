"""Tests for the decoupled scoring stage: GPU lease is released before the
CPU-bound scoring runs, and scoring can be delegated to a separate process
pool so it neither holds a GPU slot nor serializes on the yt read lock."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from grteclyn_wrapper.core.config import resolve_example
from grteclyn_wrapper.core.episode import create_episode
from grteclyn_wrapper.core.evaluation import (
    CpuGateResult,
    Evaluation,
    ScoringRequest,
    _run_gpu_session,
    _score_episode_phase,
)
from grteclyn_wrapper.search.gpu_pool import GpuPool, GpuPoolConfig
from grteclyn_wrapper.search.grtresna_evaluation_gates import (
    GRTresnaPreEvolutionGateConfig,
)
from grteclyn_wrapper.search.scoring_pool import ScoringPool


@pytest.fixture
def cpu_result(tmp_path: Path):
    episode = create_episode(tmp_path, name="scoring_episode")
    gridinit = tmp_path / "initial_data.gridinit"
    gridinit.write_bytes(b"x")
    return CpuGateResult(
        episode=episode,
        gte_overrides={"recipe_initial_data_file": str(gridinit)},
        gridinit_path=gridinit,
        gate_config=GRTresnaPreEvolutionGateConfig(),
    )


def _run(cpu_result, *, gpu_pool=None, scoring_pool=None):
    return _run_gpu_session(
        cpu_result=cpu_result,
        example=resolve_example("RadialRecipe"),
        template=resolve_example("RadialRecipe").template,
        executable=MagicMock(),
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
        grtresna_postload_gate=False,
        postload_gate_config=None,
        cleanup_plotfiles=False,
        gpu_pool=gpu_pool,
        scoring_pool=scoring_pool,
    )


def test_lease_released_before_scoring(cpu_result) -> None:
    """The GPU slot must be free by the time scoring runs."""
    pool = GpuPool(GpuPoolConfig(gpu_ids=[0], slots_per_gpu=1))
    seen_active = []

    class RecordingScorer:
        def score(self, request: ScoringRequest) -> Evaluation:
            # Lease should already be released when scoring is invoked.
            seen_active.append(pool.active_leases)
            return Evaluation(
                score=3.0, components={}, notes=[],
                episode_path=request.episode_path, exit_code=request.exit_code,
                preflight_rejected=False, reason=None, metrics={},
            )

    with patch(
        "grteclyn_wrapper.core.evaluation.reject_postload_gate", return_value=None,
    ), patch(
        "grteclyn_wrapper.core.evaluation.run_episode",
        return_value=MagicMock(returncode=0),
    ), patch("grteclyn_wrapper.core.evaluation.write_params"):
        result = _run(cpu_result, gpu_pool=pool, scoring_pool=RecordingScorer())

    assert seen_active == [0]  # no lease held during scoring
    assert pool.active_leases == 0  # fully released afterwards
    assert result.score == 3.0


def test_scoring_delegated_to_pool_not_inline(cpu_result) -> None:
    scorer = MagicMock()
    scorer.score.return_value = Evaluation(
        score=1.0, components={}, notes=[], episode_path=None, exit_code=0,
        preflight_rejected=False, reason=None, metrics={},
    )
    with patch(
        "grteclyn_wrapper.core.evaluation.reject_postload_gate", return_value=None,
    ), patch(
        "grteclyn_wrapper.core.evaluation.run_episode",
        return_value=MagicMock(returncode=0),
    ), patch("grteclyn_wrapper.core.evaluation.write_params"), patch(
        "grteclyn_wrapper.core.evaluation._score_episode_phase",
    ) as inline_score:
        _run(cpu_result, scoring_pool=scorer)

    inline_score.assert_not_called()
    scorer.score.assert_called_once()
    (request,), _ = scorer.score.call_args
    assert isinstance(request, ScoringRequest)
    assert request.episode_path == str(cpu_result.episode.path)


def test_inline_scoring_when_no_pool(cpu_result) -> None:
    with patch(
        "grteclyn_wrapper.core.evaluation.reject_postload_gate", return_value=None,
    ), patch(
        "grteclyn_wrapper.core.evaluation.run_episode",
        return_value=MagicMock(returncode=0),
    ), patch("grteclyn_wrapper.core.evaluation.write_params"), patch(
        "grteclyn_wrapper.core.evaluation._score_episode_phase",
        return_value=Evaluation(
            score=2.0, components={}, notes=[], episode_path=None, exit_code=0,
            preflight_rejected=False, reason=None, metrics={},
        ),
    ) as inline_score:
        result = _run(cpu_result)

    inline_score.assert_called_once()
    assert result.score == 2.0


def test_score_episode_phase_is_pure_function(tmp_path: Path) -> None:
    """Scoring only needs the on-disk episode dir + the request payload."""
    episode = create_episode(tmp_path, name="pure_episode")
    request = ScoringRequest(
        episode_path=str(episode.path),
        exit_code=0,
        dry_run=False,
        target_stop_time=1.0,
        score_weights=None,
        objective_mode="weighted",
        ftl_L=None,
        evolving_enabled=False,
        gte_overrides={},
        cleanup_plotfiles=False,
    )
    metrics = MagicMock(__dataclass_fields__={})
    with patch(
        "grteclyn_wrapper.core.evaluation.read_episode_metrics", return_value=metrics,
    ), patch(
        "grteclyn_wrapper.core.evaluation.score_episode",
        return_value=MagicMock(total=7.5, components={"a": 1.0}, notes=["n"]),
    ), patch(
        "grteclyn_wrapper.core.evaluation.domain_half_width_for_episode",
        return_value=1.0,
    ), patch(
        "grteclyn_wrapper.core.evaluation.dataclass_to_dict",
        side_effect=lambda obj: {"total": 7.5} if hasattr(obj, "total") else {},
    ):
        result = _score_episode_phase(request)

    assert result.score == 7.5
    assert result.components == {"a": 1.0}
    assert result.episode_path == str(episode.path)
    assert episode.score_path.is_file()


def test_scoring_pool_construction_and_bounds() -> None:
    pool = ScoringPool(max_workers=3)
    try:
        assert pool.max_workers == 3
    finally:
        pool.shutdown()

    with pytest.raises(ValueError):
        ScoringPool(max_workers=0)


def test_scoring_pool_scores_in_subprocess(tmp_path: Path) -> None:
    """End-to-end: the request is pickled to a spawned worker, scored, and an
    Evaluation is returned with score.json written to disk."""
    episode = create_episode(tmp_path, name="subproc_episode")
    request = ScoringRequest(
        episode_path=str(episode.path),
        exit_code=0,
        dry_run=False,
        target_stop_time=1.0,
        score_weights=None,
        objective_mode="weighted",
        ftl_L=None,
        evolving_enabled=False,
        gte_overrides={},
        cleanup_plotfiles=False,
    )
    pool = ScoringPool(max_workers=1)
    try:
        result = pool.score(request)
    finally:
        pool.shutdown()

    assert isinstance(result, Evaluation)
    assert result.episode_path == str(episode.path)
    assert episode.score_path.is_file()
