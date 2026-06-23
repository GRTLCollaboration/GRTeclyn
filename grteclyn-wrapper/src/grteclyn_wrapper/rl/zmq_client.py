from __future__ import annotations

import struct
from typing import Sequence

import numpy as np

try:
    import zmq
except ImportError:  # pragma: no cover - optional rl extra
    zmq = None  # type: ignore[assignment]


class ZmqObservationSource:
    def __init__(self, port: int = 5555, timeout_ms: int = 30_000) -> None:
        if zmq is None:
            raise ImportError("pyzmq is required for RL support")
        self._ctx = zmq.Context.instance()
        self._socket = self._ctx.socket(zmq.REQ)
        self._socket.setsockopt(zmq.RCVTIMEO, timeout_ms)
        self._socket.setsockopt(zmq.SNDTIMEO, timeout_ms)
        self._socket.connect(f"tcp://127.0.0.1:{port}")
        self._pending_hello = True

    def recv_obs(self) -> np.ndarray:
        if self._pending_hello:
            self._socket.send(b"\x00")
            self._pending_hello = False
        msg = self._socket.recv()
        return np.frombuffer(msg, dtype=np.float64).copy()

    def send_action(self, action: np.ndarray) -> None:
        payload = np.asarray(action, dtype=np.float64).tobytes()
        self._socket.send(payload)

    def close(self) -> None:
        self._socket.close(linger=0)


def pack_obs(obs: Sequence[float]) -> bytes:
    return struct.pack(f"{len(obs)}d", *obs)


def unpack_obs(data: bytes) -> np.ndarray:
    count = len(data) // 8
    return np.array(struct.unpack(f"{count}d", data), dtype=np.float64)
