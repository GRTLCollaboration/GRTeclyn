"""Tests for EvalPipeline scheduling."""

from __future__ import annotations

import threading
import time

from grteclyn_wrapper.core.evaluation import Evaluation
from grteclyn_wrapper.search.eval_pipeline import EvalJob, EvalPipeline
from grteclyn_wrapper.search.gpu_pool import CpuAdmissionController, GpuPool, GpuPoolConfig


def _evaluation(score: float = 1.0) -> Evaluation:
    return Evaluation(
        score=score,
        components={},
        notes=[],
        episode_path=None,
        exit_code=0,
        preflight_rejected=False,
        reason=None,
        metrics={},
    )


def test_workers_sleep_when_queue_empty_during_in_flight() -> None:
    """Idle workers must not busy-spin while GPU jobs are in flight."""
    pool = GpuPool(GpuPoolConfig(gpu_ids=[0], slots_per_gpu=1))
    admission = CpuAdmissionController(1)
    counter = [0]
    counter_lock = threading.Lock()
    started = threading.Event()
    release = threading.Event()

    def evaluate(job: EvalJob[int]) -> Evaluation:
        started.set()
        release.wait(timeout=2)
        return _evaluation()

    pipeline = EvalPipeline(
        gpu_pool=pool,
        cpu_admission=admission,
        evaluate_fn=evaluate,
        on_complete=lambda _job, _res: None,
        counter_lock=counter_lock,
        eval_counter=counter,
        max_in_flight=2,
    )
    pipeline.submit(1)
    assert started.wait(timeout=2)
    time.sleep(0.15)
    # Second worker should sleep while queue is empty and first job is in flight.
    pipeline.submit(2)
    release.set()
    pipeline.wait_until_idle()
    pipeline.shutdown()


def test_pipeline_fills_all_slots() -> None:
    pool = GpuPool(GpuPoolConfig(gpu_ids=[0, 1], slots_per_gpu=2))
    admission = CpuAdmissionController(2)
    counter = [0]
    counter_lock = threading.Lock()
    peak_in_flight = 0
    lock = threading.Lock()
    completed: list[int] = []

    def evaluate(job: EvalJob[int]) -> Evaluation:
        with lock:
            nonlocal peak_in_flight
            peak_in_flight = max(peak_in_flight, pipeline._in_flight)  # noqa: SLF001
        time.sleep(0.02)
        return _evaluation(float(job.payload))

    def on_complete(job: EvalJob[int], result: Evaluation) -> None:
        completed.append(job.eval_id)

    pipeline = EvalPipeline(
        gpu_pool=pool,
        cpu_admission=admission,
        evaluate_fn=evaluate,
        on_complete=on_complete,
        counter_lock=counter_lock,
        eval_counter=counter,
        max_in_flight=4,
    )
    pipeline.submit_many(list(range(8)))
    pipeline.wait_until_idle()
    pipeline.shutdown()
    assert peak_in_flight >= 4
    assert len(completed) == 8
    assert len(set(completed)) == 8


def test_cross_batch_backfill_unique_ids() -> None:
    pool = GpuPool(GpuPoolConfig(gpu_ids=[0], slots_per_gpu=2))
    admission = CpuAdmissionController(1)
    counter = [0]
    counter_lock = threading.Lock()
    ids: list[int] = []

    pipeline = EvalPipeline(
        gpu_pool=pool,
        cpu_admission=admission,
        evaluate_fn=lambda job: _evaluation(),
        on_complete=lambda job, _res: ids.append(job.eval_id),
        counter_lock=counter_lock,
        eval_counter=counter,
        max_in_flight=2,
    )
    pipeline.submit_many([1, 2, 3, 4, 5])
    pipeline.wait_until_idle()
    pipeline.shutdown()
    assert len(ids) == 5
    assert len(set(ids)) == 5
