"""Cross-batch evaluation pipeline with CPU admission and GPU slot leasing."""

from __future__ import annotations

import json
import threading
from collections import deque
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Generic, TypeVar

from ..core.evaluation import Evaluation
from .gpu_pool import CpuAdmissionController, GpuPool

T = TypeVar("T")


@dataclass(frozen=True)
class EvalJob(Generic[T]):
    payload: T
    eval_id: int


@dataclass
class PipelineStats:
    submitted: int = 0
    completed: int = 0
    rejected_cpu: int = 0
    rejected_gpu: int = 0


class EvalPipeline(Generic[T]):
    """Bounded sliding-window evaluator with CPU/GPU backpressure."""

    def __init__(
        self,
        *,
        gpu_pool: GpuPool,
        cpu_admission: CpuAdmissionController,
        evaluate_fn: Callable[[EvalJob[T]], Evaluation],
        on_complete: Callable[[EvalJob[T], Evaluation], None],
        counter_lock: threading.Lock,
        eval_counter: list[int],
        max_in_flight: int | None = None,
        metadata_path: Path | None = None,
        persist_counter_interval: int = 10,
    ) -> None:
        self._gpu_pool = gpu_pool
        self._cpu_admission = cpu_admission
        self._evaluate_fn = evaluate_fn
        self._on_complete = on_complete
        self._counter_lock = counter_lock
        self._eval_counter = eval_counter
        self._max_in_flight = max_in_flight or (
            gpu_pool.total_slots + cpu_admission.max_concurrent
        )
        self._metadata_path = metadata_path
        self._persist_counter_interval = persist_counter_interval
        self._queue: deque[T] = deque()
        self._in_flight = 0
        self._lock = threading.Lock()
        self._condition = threading.Condition(self._lock)
        self._workers: list[threading.Thread] = []
        self._shutdown = False
        self._stats = PipelineStats()
        self._in_flight_eval_ids: set[int] = set()

    @property
    def stats(self) -> PipelineStats:
        with self._lock:
            return PipelineStats(
                submitted=self._stats.submitted,
                completed=self._stats.completed,
                rejected_cpu=self._stats.rejected_cpu,
                rejected_gpu=self._stats.rejected_gpu,
            )

    @property
    def in_flight_eval_ids(self) -> set[int]:
        with self._lock:
            return set(self._in_flight_eval_ids)

    @property
    def cpu_admission_active(self) -> int:
        return self._cpu_admission.active

    @property
    def max_in_flight(self) -> int:
        return self._max_in_flight

    def submit(self, payload: T) -> None:
        with self._condition:
            self._queue.append(payload)
            self._stats.submitted += 1
            self._condition.notify_all()
        self._ensure_workers()

    def submit_many(self, payloads: list[T]) -> None:
        for payload in payloads:
            self.submit(payload)

    def wait_until_idle(self) -> None:
        with self._condition:
            while self._in_flight > 0 or self._queue:
                self._condition.wait()

    def shutdown(self) -> None:
        with self._condition:
            self._shutdown = True
            self._condition.notify_all()
        for worker in self._workers:
            worker.join(timeout=0.1)

    def _ensure_workers(self) -> None:
        with self._lock:
            needed = self._max_in_flight - len(self._workers)
            for idx in range(max(0, needed)):
                worker = threading.Thread(
                    target=self._worker_loop,
                    name=f"eval-pipeline-{idx}",
                    daemon=True,
                )
                self._workers.append(worker)
                worker.start()

    def _next_eval_id(self) -> int:
        with self._counter_lock:
            self._eval_counter[0] += 1
            eval_id = self._eval_counter[0]
            if (
                self._metadata_path is not None
                and eval_id % self._persist_counter_interval == 0
            ):
                self._persist_last_eval_counter(eval_id)
        return eval_id

    def _persist_last_eval_counter(self, eval_id: int) -> None:
        if self._metadata_path is None or not self._metadata_path.is_file():
            return
        try:
            metadata = json.loads(self._metadata_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            metadata = {}
        metadata["last_eval_counter"] = eval_id
        self._metadata_path.write_text(
            json.dumps(metadata, indent=2) + "\n",
            encoding="utf-8",
        )

    def _worker_loop(self) -> None:
        while True:
            with self._condition:
                while not self._shutdown and not self._queue:
                    self._condition.wait(timeout=0.1)
                if self._shutdown and not self._queue:
                    if self._in_flight == 0:
                        return
                if not self._queue or self._in_flight >= self._max_in_flight:
                    continue
                payload = self._queue.popleft()
                self._in_flight += 1

            eval_id = self._next_eval_id()
            job = EvalJob(payload=payload, eval_id=eval_id)
            with self._lock:
                self._in_flight_eval_ids.add(eval_id)

            try:
                result = self._evaluate_fn(job)
            except Exception as exc:  # noqa: BLE001 - keep pipeline alive
                result = Evaluation(
                    score=0.0,
                    components={},
                    notes=[repr(exc)],
                    episode_path=None,
                    exit_code=1,
                    preflight_rejected=False,
                    reason=repr(exc),
                    metrics={},
                )

            try:
                self._on_complete(job, result)
            finally:
                with self._condition:
                    self._in_flight -= 1
                    self._stats.completed += 1
                    self._in_flight_eval_ids.discard(eval_id)
                    self._condition.notify_all()


def scan_partial_eval_dirs(
    search_dir: Path,
    *,
    trajectory_eval_ids: set[int],
    remove_partial: bool = False,
) -> list[dict[str, Any]]:
    """Detect partial eval directories from a crashed pipeline."""
    records: list[dict[str, Any]] = []
    for child in sorted(search_dir.glob("eval_*")):
        if not child.is_dir():
            continue
        try:
            eval_id = int(child.name.split("_", 1)[1])
        except (IndexError, ValueError):
            continue
        score_path = child / "score.json"
        if score_path.is_file() or eval_id in trajectory_eval_ids:
            continue
        record = {
            "eval": eval_id,
            "status": "pipeline_interrupted",
            "episode": str(child),
            "reason": "partial_eval_without_score",
        }
        records.append(record)
        if remove_partial:
            import shutil

            from ..core.scratch import purge_plotfile_scratch

            # A crashed eval left plotfiles on node-local scratch that no
            # consumer is coming back for.  Dropping the eval dir without them
            # would leak the disk a resumed campaign is about to need.
            purge_plotfile_scratch(child, force=True)
            shutil.rmtree(child, ignore_errors=True)
    return records


def resume_eval_counter(
    *,
    metadata_path: Path | None,
    trajectory_eval_ids: set[int],
) -> int:
    meta_counter = 0
    if metadata_path is not None and metadata_path.is_file():
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            meta_counter = int(metadata.get("last_eval_counter", 0))
        except (OSError, json.JSONDecodeError, TypeError, ValueError):
            meta_counter = 0
    traj_counter = max(trajectory_eval_ids) if trajectory_eval_ids else 0
    return max(meta_counter, traj_counter)
