#!/usr/bin/env python3
"""Dummy RL agent for Gate 0B — sends zero physical actions over ZMQ."""

from __future__ import annotations

import argparse
import struct

try:
    import zmq
except ImportError as exc:  # pragma: no cover
    raise SystemExit("pyzmq required: uv sync --extra rl") from exc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=5555)
    parser.add_argument("--steps", type=int, default=100)
    args = parser.parse_args()

    ctx = zmq.Context.instance()
    sock = ctx.socket(zmq.REQ)
    sock.connect(f"tcp://127.0.0.1:{args.port}")

    zero_action = struct.pack("5d", 0.0, 0.0, 0.0, 1.0, 0.75)
    for step in range(args.steps):
        obs = sock.recv()
        print(f"step={step} obs_dim={len(obs)//8}")
        sock.send(zero_action)


if __name__ == "__main__":
    main()
