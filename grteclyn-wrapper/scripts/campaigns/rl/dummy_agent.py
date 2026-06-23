#!/usr/bin/env python3
"""Dummy RL agent for Gate 0B — sends neutral (no-op) actions over ZMQ.

Actions are raw-normalized in [-1, 1]; C++ maps them to physical bounds.  A
neutral action keeps the simulation on its IVP trajectory:
  * per lump: amp=-1 (=> 0 amplitude after C++ mapping), freq=-1, phase=0
  * gauge: lapse=0, shift=0 (=> default lapse_advec=1.0, shift_Gamma=0.75)

The observation layout (6 global + 8 per lump) lets the agent infer the lump
count, hence the expected action dimension 3*N + 2.
"""

from __future__ import annotations

import argparse
import struct

try:
    import zmq
except ImportError as exc:  # pragma: no cover
    raise SystemExit("pyzmq required: uv sync --extra rl") from exc

GLOBAL_OBS_DIM = 6
LUMP_OBS_STRIDE = 8


def _num_lumps_from_obs(num_doubles: int) -> int:
    n = (num_doubles - GLOBAL_OBS_DIM) // LUMP_OBS_STRIDE
    return max(1, n)


def _neutral_action(num_lumps: int) -> bytes:
    values: list[float] = []
    for _ in range(num_lumps):
        values += [-1.0, -1.0, 0.0]  # amp->0, freq->0, phase 0
    values += [0.0, 0.0]  # neutral gauge
    return struct.pack(f"{len(values)}d", *values)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=5555)
    parser.add_argument("--steps", type=int, default=100)
    args = parser.parse_args()

    ctx = zmq.Context.instance()
    sock = ctx.socket(zmq.REQ)
    sock.connect(f"tcp://127.0.0.1:{args.port}")

    for step in range(args.steps):
        obs = sock.recv()
        num_doubles = len(obs) // 8
        num_lumps = _num_lumps_from_obs(num_doubles)
        print(f"step={step} obs_dim={num_doubles} num_lumps={num_lumps}")
        sock.send(_neutral_action(num_lumps))


if __name__ == "__main__":
    main()
