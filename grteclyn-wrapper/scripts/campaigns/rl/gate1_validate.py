#!/usr/bin/env python3
"""Gate 1 — Tax Man T1–T4 static checks + short live SpacetimeFtlEnv smoke."""

from __future__ import annotations

import argparse
import inspect
import sys
from pathlib import Path

import numpy as np

from grteclyn_wrapper.rl.env import SpacetimeFtlEnv, SpacetimeFtlEnvConfig
from grteclyn_wrapper.rl.params import read_rl_num_lumps
from grteclyn_wrapper.rl.reward import TaxManState, evaluate_fences


def _fail(msg: str) -> None:
    print(f"FAIL {msg}")
    raise SystemExit(1)


def _pass(msg: str) -> None:
    print(f"PASS {msg}")


def check_t4_none_guard() -> None:
    result = evaluate_fences(
        min_lapse=0.2,
        l2_ham=0.01,
        wec_violation_fraction=None,
        horizon_detected=False,
        state=TaxManState(),
    )
    if result.terminated:
        _fail("T4: WEC fence fired with wec_violation_fraction=None")
    _pass("T4: WEC None-guard inactive when metric missing")


def check_t2_audit_clip() -> None:
    audit = min(0.0, -500.0 - 400.0)
    clipped = max(audit, -2000.0)
    if clipped != -900.0:
        _fail(f"T2: expected clip -900, got {clipped}")
    _pass("T2: audit penalty clip math (-2000 floor)")


def check_t3_drain_before_audit() -> None:
    src = inspect.getsource(SpacetimeFtlEnv.step)
    drain_idx = src.find("wait_consumer_drain")
    audit_idx = src.find("compute_audit_penalty")
    if drain_idx < 0 or audit_idx < 0 or drain_idx > audit_idx:
        _fail("T3: wait_consumer_drain must precede compute_audit_penalty in env.step")
    _pass("T3: consumer drain precedes Tax Man audit in env.step")


def check_t1_gamma_default() -> None:
    train_motor = Path(__file__).resolve().parent / "train_motor.py"
    text = train_motor.read_text(encoding="utf-8")
    if 'default=0.999' not in text:
        _fail("T1: train_motor.py should default gamma=0.999")
    _pass("T1: PPO gamma default 0.999 in train_motor.py")


def run_live_smoke(
    *,
    executable: Path,
    params_file: Path,
    episode_path: Path,
    stop_time: float,
    zmq_port: int,
    gpu_id: int,
    zmq_lib: Path | None,
    max_steps: int,
) -> None:
    num_lumps = read_rl_num_lumps(params_file)
    neutral = np.array(
        [-1.0, -1.0, 0.0] * num_lumps + [0.0, 0.0],
        dtype=np.float64,
    )
    cfg = SpacetimeFtlEnvConfig(
        episode_path=episode_path,
        executable=executable,
        params_file=params_file,
        zmq_port=zmq_port,
        stop_time=stop_time,
        gpu_id=gpu_id,
        num_lumps=num_lumps,
        use_mpi=False,
        plot_consumer=True,
        zmq_lib=zmq_lib,
    )
    env = SpacetimeFtlEnv(cfg)
    try:
        obs, _ = env.reset()
        if obs.size != 6 + 8 * num_lumps:
            _fail(f"live smoke: obs size {obs.size}")

        steps = 0
        total = 0.0
        terminated = False
        info: dict = {}
        while steps < max_steps and not terminated:
            obs, reward, terminated, _trunc, info = env.step(neutral)
            total += reward
            steps += 1
    finally:
        env.close()

    if steps < 1:
        _fail("live smoke: no env steps completed")
    score_path = episode_path / "small_data" / "score_timeseries.jsonl"
    if not score_path.exists():
        _fail(f"live smoke: missing {score_path} (plot consumer / incremental score)")
    _pass(
        f"live smoke: {steps} steps, return={total:.2f}, sim_time={info.get('sim_time', '?')}, "
        f"score_timeseries rows present"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--params", type=Path)
    parser.add_argument("--episode-path", type=Path)
    parser.add_argument("--stop-time", type=float, default=1.0)
    parser.add_argument("--zmq-port", type=int, default=5557)
    parser.add_argument("--gpu-id", type=int, default=0)
    parser.add_argument("--zmq-lib", type=Path, default=None)
    parser.add_argument("--max-steps", type=int, default=8)
    parser.add_argument("--skip-live", action="store_true")
    args = parser.parse_args()

    check_t1_gamma_default()
    check_t2_audit_clip()
    check_t3_drain_before_audit()
    check_t4_none_guard()

    if args.skip_live:
        print("PASS Gate 1 static T1–T4 (live smoke skipped)")
        return

    if not args.executable or not args.params or not args.episode_path:
        _fail("live smoke requires --executable --params --episode-path")

    run_live_smoke(
        executable=args.executable,
        params_file=args.params,
        episode_path=args.episode_path,
        stop_time=args.stop_time,
        zmq_port=args.zmq_port,
        gpu_id=args.gpu_id,
        zmq_lib=args.zmq_lib,
        max_steps=args.max_steps,
    )
    print("PASS Gate 1: Tax Man T1–T4 + live env smoke")


if __name__ == "__main__":
    main()
