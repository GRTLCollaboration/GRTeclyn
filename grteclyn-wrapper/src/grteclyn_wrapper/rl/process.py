from __future__ import annotations

import atexit
import os
import signal
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass
class SubprocessEpisodeLauncher:
    executable: Path
    params_file: Path
    work_dir: Path
    gpu_id: int | None = None
    mpi_ranks: int = 1

    def __post_init__(self) -> None:
        self._process: subprocess.Popen[str] | None = None
        atexit.register(self.stop)

    def start(self) -> None:
        self.stop()
        env = os.environ.copy()
        if self.gpu_id is not None:
            env["CUDA_VISIBLE_DEVICES"] = str(self.gpu_id)
        cmd = ["mpirun", "-n", str(self.mpi_ranks), str(self.executable), str(self.params_file)]
        self._process = subprocess.Popen(
            cmd,
            cwd=str(self.work_dir),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

    def stop(self) -> None:
        if self._process is None:
            return
        if self._process.poll() is None:
            self._process.send_signal(signal.SIGKILL)
            self._process.wait(timeout=5)
        self._process = None
