"""Out-of-process pool for the CPU-bound episode scoring stage.

The GPU evolution binary is fast to launch but the post-run scoring (yt reads
plus geodesic integration) is CPU-heavy and, in-process, is serialized by the
process-wide ``PLOTFILE_READ_LOCK`` because yt is not thread-safe.  Running
each episode's scoring in its own worker process:

* releases the GPU slot the moment the evolution binary exits (the caller only
  holds the lease around the binary, not the scoring), so the GPUs stay fed;
* gives each scorer its own yt state, so scoring across evals runs truly in
  parallel on the idle CPU cores instead of queueing on one lock.

Concurrency is bounded by ``max_workers`` so a burst of finished evolutions
cannot flood the node's CPUs.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

from ..core.evaluation import Evaluation, ScoringRequest, _score_episode_phase


def scoring_workers_from_env(default: int) -> int:
    """Resolve the scoring-process count, overridable via ``SCORING_WORKERS``."""
    raw = os.environ.get("SCORING_WORKERS", "").strip()
    if not raw:
        return max(1, default)
    try:
        value = int(raw)
    except ValueError:
        return max(1, default)
    return max(1, value)


def _score_worker(request: ScoringRequest) -> Evaluation:
    """Top-level entry point executed inside a scoring subprocess."""
    return _score_episode_phase(request)


class ScoringPool:
    """Bounded process pool that scores episodes off the GPU-lease path."""

    def __init__(self, max_workers: int) -> None:
        if max_workers < 1:
            raise ValueError("max_workers must be >= 1")
        self._max_workers = max_workers
        # "spawn" avoids fork-after-threads hazards (the driver runs a threaded
        # evaluation pipeline) and gives each worker a clean yt import.
        self._executor = ProcessPoolExecutor(
            max_workers=max_workers,
            mp_context=get_context("spawn"),
        )

    @property
    def max_workers(self) -> int:
        return self._max_workers

    def score(self, request: ScoringRequest) -> Evaluation:
        """Score one episode in a worker process, blocking for the result."""
        future = self._executor.submit(_score_worker, request)
        return future.result()

    def shutdown(self, *, wait: bool = True) -> None:
        self._executor.shutdown(wait=wait)

    def __enter__(self) -> "ScoringPool":
        return self

    def __exit__(self, *exc: object) -> None:
        self.shutdown()
