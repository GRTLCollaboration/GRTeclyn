from __future__ import annotations

from typing import Protocol

import numpy as np


class ObservationSource(Protocol):
    def recv_obs(self) -> np.ndarray: ...


class ActionSink(Protocol):
    def send_action(self, action: np.ndarray) -> None: ...


class EpisodeLauncher(Protocol):
    def start(self) -> None: ...

    def stop(self) -> None: ...
