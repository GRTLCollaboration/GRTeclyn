"""Unit tests for GpuPool and CPU admission."""

from __future__ import annotations

import threading
import time

import pytest

from grteclyn_wrapper.search.gpu_pool import (
    CpuAdmissionController,
    GpuPool,
    GpuPoolConfig,
    cluster_cpu_budget,
    max_concurrent_grtresna,
)


def test_multi_slot_same_gpu() -> None:
    pool = GpuPool(GpuPoolConfig(gpu_ids=[0], slots_per_gpu=5))
    active: list[int] = []
    lock = threading.Lock()

    def worker() -> None:
        with pool.lease():
            with lock:
                active.append(pool.active_leases)
            time.sleep(0.05)

    threads = [threading.Thread(target=worker) for _ in range(5)]
    for thread in threads:
        thread.start()
    time.sleep(0.02)
    with lock:
        peak = max(active) if active else 0
    for thread in threads:
        thread.join()
    assert peak == 5


def test_total_slots_8_gpus_5_slots() -> None:
    pool = GpuPool(GpuPoolConfig(gpu_ids=list(range(8)), slots_per_gpu=5))
    assert pool.total_slots == 40


def test_lease_blocks_when_pool_full() -> None:
    pool = GpuPool(GpuPoolConfig(gpu_ids=[0], slots_per_gpu=1))
    entered = threading.Event()
    release = threading.Event()

    def holder() -> None:
        with pool.lease():
            entered.set()
            release.wait(timeout=2)

    t1 = threading.Thread(target=holder)
    t1.start()
    assert entered.wait(timeout=2)

    acquired = threading.Event()

    def waiter() -> None:
        with pool.lease():
            acquired.set()

    t2 = threading.Thread(target=waiter)
    t2.start()
    time.sleep(0.05)
    assert not acquired.is_set()
    release.set()
    assert acquired.wait(timeout=2)
    t1.join(timeout=2)
    t2.join(timeout=2)


def test_grtresna_admission_blocks_excess() -> None:
    admission = CpuAdmissionController(3)
    peak = 0
    lock = threading.Lock()

    def worker() -> None:
        nonlocal peak
        with admission.admit():
            with lock:
                peak = max(peak, admission.active)
            time.sleep(0.05)

    threads = [threading.Thread(target=worker) for _ in range(10)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert peak == 3


def test_release_on_exception() -> None:
    pool = GpuPool(GpuPoolConfig(gpu_ids=[0], slots_per_gpu=1))
    with pytest.raises(RuntimeError):
        with pool.lease():
            raise RuntimeError("boom")
    with pool.lease():
        pass


def test_cluster_cpu_budget_and_max_grtresna() -> None:
    assert cluster_cpu_budget(192, cluster_cpu_fraction=0.30) == 57
    assert max_concurrent_grtresna(58, mpi_ranks=8, reserve_cores=4) == 6
