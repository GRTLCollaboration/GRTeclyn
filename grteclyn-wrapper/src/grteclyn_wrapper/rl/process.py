from __future__ import annotations

import atexit
import os
import signal
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class SubprocessEpisodeLauncher:
    executable: Path
    params_file: Path
    work_dir: Path
    gpu_id: int | None = None
    mpi_ranks: int = 1
    use_mpi: bool = False
    plot_consumer: bool = False
    plot_consumer_frames: bool = False
    plot_consumer_jobs: int = 1
    objective_mode: str = "general_ftl"
    target_stop_time: float | None = None
    zmq_lib: Path | None = None
    _process: subprocess.Popen[str] | None = field(default=None, init=False, repr=False)
    _consumer: subprocess.Popen[str] | None = field(default=None, init=False, repr=False)
    _sim_log: Any = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        atexit.register(self.stop)

    def _build_env(self) -> dict[str, str]:
        env = os.environ.copy()
        if self.gpu_id is not None:
            env["CUDA_VISIBLE_DEVICES"] = str(self.gpu_id)
        ld_parts: list[str] = []
        if self.zmq_lib is not None:
            ld_parts.append(str(self.zmq_lib.resolve()))
        for part in env.get("LD_LIBRARY_PATH", "").split(":"):
            if not part:
                continue
            # Conda/grtresna libstdc++ can ABI-crash the CUDA binary during AMReX regrid.
            if ".mlspace/envs/" in part or part.endswith("/grtresna/lib"):
                continue
            ld_parts.append(part)
        if ld_parts:
            env["LD_LIBRARY_PATH"] = ":".join(dict.fromkeys(ld_parts))
        return env

    def sim_running(self) -> bool:
        return self._process is not None and self._process.poll() is None

    def sim_returncode(self) -> int | None:
        if self._process is None:
            return None
        return self._process.poll()

    def _build_command(self) -> list[str]:
        base = [str(self.executable), str(self.params_file)]
        if self.use_mpi and self.mpi_ranks > 1:
            return ["mpirun", "-n", str(self.mpi_ranks), *base]
        return base

    def start(self) -> None:
        self.stop()
        self.work_dir.mkdir(parents=True, exist_ok=True)

        if self.plot_consumer:
            from grteclyn_wrapper.core.episode import Episode
            from grteclyn_wrapper.core.runner import start_plotfile_consumer

            episode = Episode(self.work_dir)
            self._consumer = start_plotfile_consumer(
                episode,
                example_name="RadialRecipe",
                incremental_score=True,
                objective_mode=self.objective_mode,
                target_stop_time=self.target_stop_time,
                ftl_timeseries=True,
                frames=self.plot_consumer_frames,
                jobs=self.plot_consumer_jobs,
            )

        self._sim_log = open(self.work_dir / "sim_stdout.log", "a", encoding="utf-8")
        self._process = subprocess.Popen(
            self._build_command(),
            cwd=str(self.work_dir),
            env=self._build_env(),
            stdout=self._sim_log,
            stderr=subprocess.STDOUT,
            text=True,
        )

    def stop(self) -> None:
        if self._process is not None and self._process.poll() is None:
            self._process.send_signal(signal.SIGTERM)
            try:
                self._process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self._process.send_signal(signal.SIGKILL)
                self._process.wait(timeout=5)
        self._process = None

        if self._sim_log is not None:
            self._sim_log.close()
            self._sim_log = None

        if self._consumer is not None:
            from grteclyn_wrapper.core.runner import stop_process

            stop_process(self._consumer)
            self._consumer = None
