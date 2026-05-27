from __future__ import annotations

import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from .config import ExecutableConfig, ExampleConfig, REPO_ROOT, resolve_example
from .episode import Episode, update_metadata
from .plot_consumer import build_consume_command, default_radii_for_example


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
    example_name: str = "RadialRecipe",
    radii: Sequence[float] | None = None,
    n_points: int = 64,
    delete: bool = True,
    keep_last: int = 1,
    frames: bool = True,
    jobs: int = 4,
) -> subprocess.Popen[str]:
    profile = "wormhole" if example_name == "SupportedWormholeCollapse" else "radial"
    if radii is None:
        radii = default_radii_for_example(example_name)
    command = build_consume_command(
        episode,
        profile=profile,
        radii=radii,
        n_points=n_points,
        delete=delete,
        keep_last=keep_last,
        watch=True,
        jobs=jobs,
        frames=frames,
    )

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


def drain_plotfile_backlog(
    episode: Episode,
    *,
    example_name: str = "RadialRecipe",
    radii: Sequence[float] | None = None,
    n_points: int = 64,
    delete: bool = True,
    keep_last: int = 1,
    frames: bool = True,
    jobs: int = 4,
) -> RunResult:
    """Process any plotfiles left after the watch consumer stops."""
    profile = "wormhole" if example_name == "SupportedWormholeCollapse" else "radial"
    if radii is None:
        radii = default_radii_for_example(example_name)
    command = build_consume_command(
        episode,
        profile=profile,
        radii=radii,
        n_points=n_points,
        delete=delete,
        keep_last=keep_last,
        watch=False,
        jobs=jobs,
        frames=frames,
    )
    return _run_and_tee(command, episode.log_path, cwd=REPO_ROOT)


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
    example_dir = executable.example.dir
    if not executable.path.exists():
        raise FileNotFoundError(
            f"Executable not found: {executable.path}. Build {executable.example.name} first."
        )

    env = _merged_env(cuda_devices=cuda_devices, extra_env=extra_env)
    update_metadata(
        episode,
        {
            "executable": str(executable.path),
            "example": executable.example.name,
            "mpi_ranks": executable.mpi_ranks,
            "cuda_devices": cuda_devices,
        },
    )

    if check_params:
        check_command = build_command(executable, episode.params_path, "check_params=1")
        check_result = _run_and_tee(check_command, episode.log_path, cwd=example_dir, env=env)
        if check_result.returncode != 0:
            update_metadata(episode, {"check_params_exit_code": check_result.returncode})
            raise RuntimeError(f"check_params failed with exit code {check_result.returncode}")

    consumer = None
    if consume_plotfiles:
        consumer = start_plotfile_consumer(
            episode,
            example_name=executable.example.name,
            radii=consumer_radii,
            delete=consumer_delete,
            keep_last=1,
            frames=False,
        )

    try:
        command = build_command(executable, episode.params_path)
        result = _run_and_tee(command, episode.log_path, cwd=example_dir, env=env)
    finally:
        if consumer is not None:
            stop_process(consumer, timeout=2.0)
            drain_plotfile_backlog(
                episode,
                example_name=executable.example.name,
                radii=consumer_radii,
                delete=consumer_delete,
                keep_last=1,
                frames=True,
            )

    update_metadata(
        episode,
        {
            "simulation_exit_code": result.returncode,
            "simulation_elapsed_seconds": result.elapsed_seconds,
        },
    )
    return result
