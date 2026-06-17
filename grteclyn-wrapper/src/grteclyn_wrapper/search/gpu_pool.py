"""GPU slot pool and CPU admission control for pipelined search evaluation."""

from __future__ import annotations

import os
import threading
from collections import deque
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Iterator, Sequence


@dataclass(frozen=True)
class GpuPoolConfig:
    gpu_ids: Sequence[int]
    slots_per_gpu: int = 1
    memory_gb_per_job: float | None = None
    total_memory_gb: float | None = None

    def __post_init__(self) -> None:
        if self.slots_per_gpu < 1:
            raise ValueError("slots_per_gpu must be >= 1")
        if not self.gpu_ids:
            raise ValueError("gpu_ids must not be empty")


def cluster_cpu_budget(
    cpu_count: int,
    *,
    cluster_cpu_fraction: float = 0.30,
    pipeline_share: float = 1.0,
) -> int:
    return max(1, int(cpu_count * cluster_cpu_fraction * pipeline_share))


def max_concurrent_grtresna(
    cpu_budget: int,
    *,
    mpi_ranks: int = 8,
    reserve_cores: int = 4,
) -> int:
    return max(1, (cpu_budget - reserve_cores) // mpi_ranks)


def pipeline_share_from_env() -> float:
  raw = os.environ.get("PIPELINE_CPU_SHARE", "").strip()
  if not raw:
      return 1.0
  return float(raw)


def cluster_cpu_fraction_from_env(default: float = 0.30) -> float:
    raw = os.environ.get("CLUSTER_CPU_FRACTION", "").strip()
    if not raw:
        return default
    return float(raw)


def reserve_cores_for_pipeline(*, pipeline_share: float) -> int:
    """MODE=par branches use a smaller reserve to stay within the shared budget."""
    return 2 if pipeline_share < 1.0 else 4


class GpuPool:
    """Slot-based GPU lease pool with blocking backpressure."""

    def __init__(self, config: GpuPoolConfig) -> None:
        self._config = config
        self._slots: deque[str] = deque()
        for gpu_id in config.gpu_ids:
            device = str(gpu_id)
            for _ in range(config.slots_per_gpu):
                self._slots.append(device)
        self._condition = threading.Condition()
        self._active = 0

    @property
    def total_slots(self) -> int:
        return len(self._config.gpu_ids) * self._config.slots_per_gpu

    @property
    def slots_per_gpu(self) -> int:
        return self._config.slots_per_gpu

    @property
    def gpu_ids(self) -> Sequence[int]:
        return self._config.gpu_ids

    @property
    def active_leases(self) -> int:
        with self._condition:
            return self._active

    @contextmanager
    def lease(self) -> Iterator[str]:
        with self._condition:
            while not self._slots:
                self._condition.wait()
            device = self._slots.popleft()
            self._active += 1
        try:
            yield device
        finally:
            with self._condition:
                self._slots.append(device)
                self._active -= 1
                self._condition.notify()


class CpuAdmissionController:
    """Semaphore limiting concurrent GRTresna CPU-phase work."""

    def __init__(self, max_concurrent: int) -> None:
        if max_concurrent < 1:
            raise ValueError("max_concurrent must be >= 1")
        self._max_concurrent = max_concurrent
        self._semaphore = threading.Semaphore(max_concurrent)
        self._active = 0
        self._lock = threading.Lock()

    @property
    def max_concurrent(self) -> int:
        return self._max_concurrent

    @property
    def active(self) -> int:
        with self._lock:
            return self._active

    @contextmanager
    def admit(self) -> Iterator[None]:
        self._semaphore.acquire()
        with self._lock:
            self._active += 1
        try:
            yield
        finally:
            with self._lock:
                self._active -= 1
            self._semaphore.release()


def build_pipeline_resources(
    gpu_ids: Sequence[int],
    *,
    slots_per_gpu: int = 1,
    cpu_count: int | None = None,
    cluster_cpu_fraction: float | None = None,
    pipeline_share: float | None = None,
    mpi_ranks: int = 8,
    reserve_cores: int | None = None,
    memory_gb_per_job: float | None = None,
) -> tuple[GpuPool, CpuAdmissionController, dict[str, int | float]]:
    """Construct GpuPool + CpuAdmissionController with startup sizing metadata."""
    cpu_total = cpu_count if cpu_count is not None else os.cpu_count() or 1
    fraction = (
        cluster_cpu_fraction
        if cluster_cpu_fraction is not None
        else cluster_cpu_fraction_from_env()
    )
    share = pipeline_share if pipeline_share is not None else pipeline_share_from_env()
    reserve = reserve_cores if reserve_cores is not None else reserve_cores_for_pipeline(
        pipeline_share=share,
    )
    budget = cluster_cpu_budget(cpu_total, cluster_cpu_fraction=fraction, pipeline_share=share)
    max_grt = max_concurrent_grtresna(budget, mpi_ranks=mpi_ranks, reserve_cores=reserve)
    gpu_pool = GpuPool(
        GpuPoolConfig(
            gpu_ids=gpu_ids,
            slots_per_gpu=slots_per_gpu,
            memory_gb_per_job=memory_gb_per_job,
        )
    )
    cpu_admission = CpuAdmissionController(max_grt)
    sizing = {
        "node_cpus": cpu_total,
        "cluster_cpu_fraction": fraction,
        "cpu_budget": budget,
        "pipeline_share": share,
        "mpi_ranks": mpi_ranks,
        "reserve_cores": reserve,
        "max_grtresna": max_grt,
        "gpu_slots": gpu_pool.total_slots,
        "slots_per_gpu": slots_per_gpu,
    }
    print(
        "[pipeline] "
        f"node_cpus={cpu_total} cluster_cpu_fraction={fraction} cpu_budget={budget} "
        f"pipeline_share={share}"
    )
    print(
        "[pipeline] "
        f"mpi_ranks={mpi_ranks} max_grtresna={max_grt} gpu_slots={gpu_pool.total_slots} "
        f"slots_per_gpu={slots_per_gpu}"
    )
    return gpu_pool, cpu_admission, sizing
