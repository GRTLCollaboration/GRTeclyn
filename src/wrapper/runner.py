from __future__ import annotations

import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from .config import ExecutableConfig, REPO_ROOT, SUPPORTED_WORMHOLE_DIR
from .episode import Episode, update_metadata


@dataclass(frozen=True)
class RunResult:
    command: list[str]
    returncode: int
    elapsed_seconds: float


def build_command(executable: ExecutableConfig, params_path: Path, *extra_args: str) -> list[str]:
    base = [str(executable.path), str(params_path), *extra_args]
    if executable.uses_mpi:
        return ["mpirun", "-n", str(executable.mpi_ranks), *base]
    return base


def _merged_env(cuda_devices: str | None = None, extra_env: Mapping[str, str] | None = None) -> dict[str, str]:
    env = dict(os.environ)
    if cuda_devices is not None:
        env["CUDA_VISIBLE_DEVICES"] = cuda_devices
    if extra_env:
        env.update({str(key): str(value) for key, value in extra_env.items()})
    return env


def _run_and_tee(
    command: Sequence[str],
    log_path: Path,
    *,
    cwd: Path,
    env: Mapping[str, str] | None = None,
) -> RunResult:
    start = time.monotonic()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n$ {' '.join(command)}\n")
        log.flush()
        process = subprocess.Popen(
            list(command),
            cwd=cwd,
            env=dict(env) if env is not None else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            log.write(line)
            log.flush()
            print(line, end="")
        returncode = process.wait()
        elapsed = time.monotonic() - start
        log.write(f"\n[wrapper] exit_code={returncode} elapsed_seconds={elapsed:.3f}\n")

    return RunResult(command=list(command), returncode=returncode, elapsed_seconds=elapsed)


def start_plotfile_consumer(
    episode: Episode,
    *,
    radii: Sequence[float] = (8.0, 16.0),
    n_points: int = 64,
    delete: bool = False,
    keep_last: int = 2,
    frames: bool = False,
) -> subprocess.Popen[str]:
    command = [
        sys.executable,
        "-m",
        "src.visualisation.process_wave.consume_plotfiles",
        "--data",
        str(episode.path),
        "--out",
        str(episode.small_data_dir),
        "--radii",
        *[str(r) for r in radii],
        "--n-points",
        str(n_points),
        "--areal-radius",
        "--watch",
    ]
    if delete:
        command.extend(["--delete", "--keep-last", str(keep_last)])
    if frames:
        command.extend(["--frames-out", str(episode.frames_dir)])

    log = episode.log_path.open("a", encoding="utf-8")
    log.write(f"\n$ {' '.join(command)}\n")
    log.flush()
    return subprocess.Popen(
        command,
        cwd=REPO_ROOT,
        stdout=log,
        stderr=subprocess.STDOUT,
        text=True,
    )


def stop_process(process: subprocess.Popen[str], timeout: float = 10.0) -> None:
    if process.poll() is not None:
        return
    process.terminate()
    try:
        process.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait(timeout=timeout)


def run_episode(
    episode: Episode,
    executable: ExecutableConfig,
    *,
    check_params: bool = True,
    cuda_devices: str | None = None,
    extra_env: Mapping[str, str] | None = None,
    consume_plotfiles: bool = False,
    consumer_radii: Sequence[float] = (8.0, 16.0),
    consumer_delete: bool = False,
) -> RunResult:
    if not executable.path.exists():
        raise FileNotFoundError(
            f"Executable not found: {executable.path}. Build SupportedWormholeCollapse first."
        )

    env = _merged_env(cuda_devices=cuda_devices, extra_env=extra_env)
    update_metadata(
        episode,
        {
            "executable": str(executable.path),
            "mpi_ranks": executable.mpi_ranks,
            "cuda_devices": cuda_devices,
        },
    )

    if check_params:
        check_command = build_command(executable, episode.params_path, "check_params=1")
        check_result = _run_and_tee(check_command, episode.log_path, cwd=SUPPORTED_WORMHOLE_DIR, env=env)
        if check_result.returncode != 0:
            update_metadata(episode, {"check_params_exit_code": check_result.returncode})
            raise RuntimeError(f"check_params failed with exit code {check_result.returncode}")

    consumer = None
    if consume_plotfiles:
        consumer = start_plotfile_consumer(
            episode,
            radii=consumer_radii,
            delete=consumer_delete,
        )

    try:
        command = build_command(executable, episode.params_path)
        result = _run_and_tee(command, episode.log_path, cwd=SUPPORTED_WORMHOLE_DIR, env=env)
    finally:
        if consumer is not None:
            stop_process(consumer)

    update_metadata(
        episode,
        {
            "simulation_exit_code": result.returncode,
            "simulation_elapsed_seconds": result.elapsed_seconds,
        },
    )
    return result
