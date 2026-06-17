"""Concurrency and resume safety tests for pipelined QD."""

from __future__ import annotations

import json
import threading
from pathlib import Path

from grteclyn_wrapper.search.eval_pipeline import resume_eval_counter, scan_partial_eval_dirs
from grteclyn_wrapper.search.qd_search.archive import Elite, QDArchive
from grteclyn_wrapper.search.qd_search.driver import run_qd_search
from grteclyn_wrapper.search.validation_tiers import DEFAULT_TIER_CONFIG, evaluate_tiers
from grteclyn_wrapper.search.qd_search.descriptors import _descriptor_details


def test_archive_lock_prevents_concurrent_insert_race() -> None:
    archive = QDArchive(bins=4)
    lock = threading.Lock()
    errors: list[Exception] = []

    def insert_worker(score: float) -> None:
        try:
            elite = Elite(
                cell=(0, 0),
                score=score,
                descriptors=(0.1, 0.2),
                params={"a": score},
                episode=None,
                tier=1,
                tier_name="test",
                objectives={},
                descriptor_details={},
            )
            with lock:
                archive.insert(elite)
        except Exception as exc:  # noqa: BLE001
            errors.append(exc)

    threads = [threading.Thread(target=insert_worker, args=(float(i),)) for i in range(8)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert not errors


def test_record_build_outside_lock_tier_computation() -> None:
    components = {"ftl": 1.0}
    metrics = {"ftl_timeseries": None}
    details = _descriptor_details(components, metrics, mode="legacy")
    assessment = evaluate_tiers(components, metrics=metrics, config=DEFAULT_TIER_CONFIG)
    assert details["x"] is not None
    assert assessment.tier_name


def test_prune_protects_in_flight_eval_ids(tmp_path: Path) -> None:
    qd_dir = tmp_path / "qd_test"
    qd_dir.mkdir()
    in_flight = {7}
    (qd_dir / "eval_000007").mkdir()
    (qd_dir / "eval_000001").mkdir()
    trajectory = [{"eval": 1, "score": 1.0, "status": "gpu_ok"}]
    from grteclyn_wrapper.search.qd_search.io import _prune_eval_dirs
    from grteclyn_wrapper.search.ftl_retention import compute_keep_eval_ids

    keep_ids = compute_keep_eval_ids(
        trajectory,
        keep_top_score=1,
        board=None,
        ftl_retention_enabled=False,
        protect_eval_ids=in_flight,
    )
    _prune_eval_dirs(qd_dir, trajectory, keep_eval_ids=keep_ids, protect_eval_ids=in_flight)
    assert (qd_dir / "eval_000007").is_dir()


def test_resume_uses_last_eval_counter_from_metadata(tmp_path: Path) -> None:
    metadata_path = tmp_path / "metadata.json"
    metadata_path.write_text(json.dumps({"last_eval_counter": 20}), encoding="utf-8")
    counter = resume_eval_counter(metadata_path=metadata_path, trajectory_eval_ids={15, 18})
    assert counter == 20


def test_scan_partial_eval_dirs_defaults_to_keep(tmp_path: Path) -> None:
    qd_dir = tmp_path / "qd"
    qd_dir.mkdir()
    partial = qd_dir / "eval_000010"
    partial.mkdir()
    records = scan_partial_eval_dirs(qd_dir, trajectory_eval_ids=set(), remove_partial=False)
    assert records and records[0]["status"] == "pipeline_interrupted"
    assert partial.is_dir()


def test_no_midrun_trajectory_compaction(tmp_path: Path, monkeypatch) -> None:
    from grteclyn_wrapper.core.evaluation import Evaluation
    from grteclyn_wrapper.search import qd_search as qd_module

    def fake_evaluate_overrides(overrides, *, out_dir, name, **_kwargs):
        episode = out_dir / name
        episode.mkdir(parents=True, exist_ok=True)
        return Evaluation(
            score=1.0,
            components={"constraint_growth": 1.0, "operational_ftl": 0.0},
            notes=[],
            episode_path=str(episode),
            exit_code=0,
            preflight_rejected=False,
            reason=None,
            metrics={},
        )

    monkeypatch.setattr(qd_module, "evaluate_overrides", fake_evaluate_overrides)

    rewrite_calls: list[int] = []

    original_open = Path.open

    def counting_open(self, *args, **kwargs):
        if self.name == "trajectory.jsonl" and args and args[0] == "w":
            rewrite_calls.append(1)
        return original_open(self, *args, **kwargs)

    monkeypatch.setattr(Path, "open", counting_open)

    run_qd_search(
        runs_dir=tmp_path,
        name="qd_no_compact",
        iterations=0,
        init_random=2,
        batch_size=2,
        bins=4,
        seed=3,
        dry_run=True,
        gpu_ids=[0],
        use_pipeline=True,
        grtresna=False,
        use_preflight=False,
        constrained=False,
        check_params=False,
        consume_plotfiles=False,
    )
    assert rewrite_calls == []


def test_qd_pipeline_dry_run_smoke(tmp_path: Path) -> None:
    archive = run_qd_search(
        runs_dir=tmp_path,
        iterations=1,
        batch_size=2,
        bins=4,
        init_random=2,
        seed=1,
        dry_run=True,
        gpu_ids=[0],
        use_pipeline=True,
        grtresna=False,
        use_preflight=False,
        constrained=False,
    )
    assert isinstance(archive, QDArchive)
