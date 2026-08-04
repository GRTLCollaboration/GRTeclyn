"""How a GRTresna solve is launched: mpirun vs. singleton.

A single rank never needs a launcher, and on a node whose PRRTE is broken
(mpirun segfaults before it launches anything) routing a one-rank solve through
mpirun turns every candidate into a spurious ``grtresna_failed`` rejection.  The
old launch test also required the executable name to lack ``.MPI.``, which no
GRTresna build on this site satisfies -- so the serial path was unreachable.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from grteclyn_wrapper.grtresna.solver import runner as solver_runner
from grteclyn_wrapper.grtresna.solver.config import GRTresnaConfig


class _FakeProc:
    """Popen stand-in: records nothing, exits cleanly, produces no output."""

    returncode = 0

    def communicate(self, timeout: int | None = None) -> tuple[str, str]:
        del timeout
        return "", ""


def _config_with_mpi_exe(tmp_path: Path, ranks: int) -> GRTresnaConfig:
    example_dir = tmp_path / "Examples" / "ScalarFieldBH"
    example_dir.mkdir(parents=True)
    exe = example_dir / "Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex"
    exe.write_text("#!/bin/sh\n", encoding="utf-8")
    exe.chmod(0o755)
    return GRTresnaConfig(grtresna_root=tmp_path, example="ScalarFieldBH", mpi_ranks=ranks)


def _capture_launch(monkeypatch, cfg: GRTresnaConfig, work_dir: Path) -> list[str]:
    captured: list[str] = []

    def fake_popen(cmd_parts, **kwargs):  # noqa: ANN001, ANN003
        del kwargs
        captured.extend(cmd_parts)
        return _FakeProc()

    monkeypatch.setattr(solver_runner.subprocess, "Popen", fake_popen)
    # The solve dies at the missing output file; the command is already built.
    with pytest.raises(FileNotFoundError):
        solver_runner.solve(cfg, work_dir=work_dir)
    return captured


def test_single_rank_runs_mpi_binary_directly(tmp_path, monkeypatch) -> None:
    """One rank launches the executable itself, even when it is an MPI build."""
    cfg = _config_with_mpi_exe(tmp_path, ranks=1)

    def explode(_cfg):  # noqa: ANN001
        raise AssertionError("a one-rank solve must not go looking for mpirun")

    monkeypatch.setattr(solver_runner, "_resolve_mpirun", explode)

    cmd = _capture_launch(monkeypatch, cfg, tmp_path / "work")

    assert cmd[0].endswith(".MPI.ex")
    assert cmd[1].endswith("params.txt")
    assert not any("mpirun" in part for part in cmd)
    assert "-np" not in cmd


def test_multi_rank_still_uses_mpirun(tmp_path, monkeypatch) -> None:
    """The launcher path is unchanged for genuine multi-rank solves."""
    cfg = _config_with_mpi_exe(tmp_path, ranks=8)
    fake_mpirun = tmp_path / "fake_mpirun"
    fake_mpirun.write_text("#!/bin/sh\n", encoding="utf-8")
    fake_mpirun.chmod(0o755)
    monkeypatch.setattr(
        solver_runner, "_resolve_mpirun", lambda _cfg: (str(fake_mpirun), {})
    )

    cmd = _capture_launch(monkeypatch, cfg, tmp_path / "work")

    assert cmd[0] == str(fake_mpirun)
    assert "-np" in cmd and cmd[cmd.index("-np") + 1] == "8"
