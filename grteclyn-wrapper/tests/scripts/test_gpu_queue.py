"""Tests for the press-and-forget work-queue runner (lib/gpu_queue.sh)."""

from __future__ import annotations

import fcntl
import os
import subprocess
from pathlib import Path

WRAPPER_ROOT = Path(__file__).resolve().parents[2]
QUEUE_SCRIPT = WRAPPER_ROOT / "scripts/campaigns/lib/gpu_queue.sh"


def _run_queue(
    qdir: Path,
    slots: list[str],
    *,
    extra_env: dict[str, str] | None = None,
    timeout: int = 60,
) -> subprocess.CompletedProcess[str]:
    env = {
        **os.environ,
        "QUEUE_POLL_SEC": "1",
        "QUEUE_IDLE_EXIT": "1",
    }
    if extra_env:
        env.update(extra_env)
    return subprocess.run(
        ["bash", str(QUEUE_SCRIPT), str(qdir), *slots],
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def _enqueue(qdir: Path, name: str, body: str) -> None:
    (qdir / "pending").mkdir(parents=True, exist_ok=True)
    (qdir / "pending" / name).write_text(body)


def test_jobs_run_failures_isolated_and_chained_jobs_picked_up(tmp_path: Path) -> None:
    qdir = tmp_path / "q"
    _enqueue(qdir, "010_ok.job", 'echo "ok on $QUEUE_SLOT"\n')
    _enqueue(qdir, "020_bad.job", "exit 7\n")
    # A "solve" job that enqueues its "evolve" follow-up — the prestage pattern.
    _enqueue(
        qdir,
        "030_chain.job",
        f'printf \'echo evolved on "$QUEUE_GPU"\\n\' > "{qdir}/pending/040_evolve.job"\n',
    )
    proc = _run_queue(qdir, ["s1", "s2"])
    assert proc.returncode == 0, proc.stderr
    done = sorted(p.name for p in (qdir / "done").iterdir())
    assert done == ["010_ok.job", "030_chain.job", "040_evolve.job"]
    failed = [p.name for p in (qdir / "failed").iterdir()]
    assert failed == ["020_bad.job"]
    assert "ok on s" in (qdir / "logs" / "010_ok.job.log").read_text()
    assert "evolved on s" in (qdir / "logs" / "040_evolve.job.log").read_text()


def test_jobs_dispatch_in_lexicographic_order(tmp_path: Path) -> None:
    qdir = tmp_path / "q"
    marker = qdir / "order.txt"
    for idx in ("030", "010", "020"):
        _enqueue(qdir, f"{idx}_j.job", f'echo {idx} >> "{marker}"\n')
    proc = _run_queue(qdir, ["s1"])  # single slot => strictly ordered
    assert proc.returncode == 0, proc.stderr
    assert marker.read_text().split() == ["010", "020", "030"]


def test_stop_sentinel_prevents_dispatch(tmp_path: Path) -> None:
    qdir = tmp_path / "q"
    _enqueue(qdir, "010_never.job", "echo should-not-run\n")
    (qdir / "STOP").touch()
    proc = _run_queue(qdir, ["s1"], extra_env={"QUEUE_IDLE_EXIT": "0"})
    assert proc.returncode == 0, proc.stderr
    assert (qdir / "pending" / "010_never.job").exists()
    assert not (qdir / "done").exists() or not list((qdir / "done").iterdir())
    assert "STOP seen" in (qdir / "queue.log").read_text()


def test_second_runner_refused_while_lock_held(tmp_path: Path) -> None:
    qdir = tmp_path / "q"
    (qdir / "pending").mkdir(parents=True)
    lock_path = qdir / "runner.lock"
    with open(lock_path, "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        proc = _run_queue(qdir, ["s1"])
        assert proc.returncode == 1
        assert "another runner already owns" in proc.stderr


def test_stale_claims_requeued_on_start(tmp_path: Path) -> None:
    qdir = tmp_path / "q"
    (qdir / "running").mkdir(parents=True)
    (qdir / "pending").mkdir(parents=True)
    stale = qdir / "running" / "010_stale.job.slots9"
    stale.write_text('echo recovered on "$QUEUE_SLOT"\n')
    proc = _run_queue(qdir, ["s1"])
    assert proc.returncode == 0, proc.stderr
    assert not stale.exists()
    assert (qdir / "done" / "010_stale.job").exists()
    assert "requeued stale claim 010_stale.job" in (qdir / "queue.log").read_text()
