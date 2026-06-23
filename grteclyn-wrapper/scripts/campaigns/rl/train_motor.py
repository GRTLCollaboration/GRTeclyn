#!/usr/bin/env python3
"""PPO trainer for Spacetime FTL RL (Phase 3)."""

from __future__ import annotations

import argparse
from pathlib import Path

from stable_baselines3 import PPO
from stable_baselines3.common.vec_env import SubprocVecEnv

from grteclyn_wrapper.rl.env import SpacetimeFtlEnv, SpacetimeFtlEnvConfig


def _make_env(port: int, gpu_id: int, episode_path: Path, executable: Path, params_file: Path):
    def _init():
        return SpacetimeFtlEnv(
            SpacetimeFtlEnvConfig(
                episode_path=episode_path,
                executable=executable,
                params_file=params_file,
                zmq_port=port,
                gpu_id=gpu_id,
            )
        )

    return _init


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--params", type=Path, required=True)
    parser.add_argument("--episode-path", type=Path, required=True)
    parser.add_argument("--n-envs", type=int, default=8)
    parser.add_argument("--base-port", type=int, default=5555)
    parser.add_argument("--timesteps", type=int, default=1_000_000)
    parser.add_argument("--gamma", type=float, default=0.999)
    parser.add_argument("--gae-lambda", type=float, default=0.95)
    args = parser.parse_args()

    env_fns = [
        _make_env(
            args.base_port + i,
            i,
            args.episode_path,
            args.executable,
            args.params,
        )
        for i in range(args.n_envs)
    ]
    env = SubprocVecEnv(env_fns)
    model = PPO(
        "MlpPolicy",
        env,
        gamma=args.gamma,
        gae_lambda=args.gae_lambda,
        verbose=1,
    )
    model.learn(total_timesteps=args.timesteps)
    model.save("spacetime_ftl_ppo")


if __name__ == "__main__":
    main()
