#!/usr/bin/env python3
"""Phase 2 baselines: zero/defensive/random gauge policies."""

from __future__ import annotations

import argparse
import math

import numpy as np

from grteclyn_wrapper.rl.reward import TaxManState, baseline_dense_reward, compute_dense_reward


def run_baseline(name: str, steps: int) -> float:
    state = TaxManState()
    total = 0.0
    for step in range(steps):
        l2_ham = 0.01 + 0.001 * step
        min_lapse = 0.25
        if name == "steer-not-break":
            total += baseline_dense_reward(l2_ham=l2_ham, min_lapse=min_lapse, max_abs_k=0.2, alpha=1.0)
        elif name == "kamikaze":
            ftl = 0.19 if step < 3 else 0.0
            total += compute_dense_reward(ftl_geo=ftl, l2_ham=l2_ham, min_lapse=min_lapse, state=state)
            if step >= 5:
                break
        else:
            total += baseline_dense_reward(l2_ham=l2_ham, min_lapse=min_lapse, max_abs_k=0.1)
    return total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=10)
    args = parser.parse_args()
    for name in ("zero", "steer-not-break", "kamikaze"):
        value = run_baseline(name, args.steps)
        print(f"{name}: return={value:.3f} frames_per_episode~={args.steps} gamma>=0.999")


if __name__ == "__main__":
    main()
