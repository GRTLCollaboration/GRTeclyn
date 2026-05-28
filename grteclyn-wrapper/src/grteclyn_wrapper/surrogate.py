"""Lightweight numpy-only RBF surrogate for surrogate-assisted search.

A full GPU evolution per evaluation is the search bottleneck.  This module
fits a cheap radial-basis-function regressor on the growing ``(theta -> S)``
archive and uses it to pre-screen CMA-ES proposals: only the most promising
(and most uncertain) candidates are evaluated on the GPU; the rest receive a
surrogate-predicted fitness so the evolution strategy still updates.

No SciPy/scikit-learn dependency: the kernel ridge solve is a small dense
linear system handled by ``numpy.linalg``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass
class RBFSurrogate:
    """Gaussian-RBF kernel-ridge regressor over a normalized search box."""

    lower: NDArray[np.float64]
    upper: NDArray[np.float64]
    ridge: float = 1.0e-6
    epsilon: float = 0.5

    _x_train: NDArray[np.float64] | None = None
    _weights: NDArray[np.float64] | None = None
    _y_mean: float = 0.0

    def _normalize(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        span = np.where((self.upper - self.lower) > 0, self.upper - self.lower, 1.0)
        return (np.asarray(x, dtype=float) - self.lower) / span

    @staticmethod
    def _pairwise_sq_dists(a: NDArray[np.float64], b: NDArray[np.float64]) -> NDArray[np.float64]:
        a2 = np.sum(a * a, axis=1)[:, None]
        b2 = np.sum(b * b, axis=1)[None, :]
        return np.maximum(a2 + b2 - 2.0 * a @ b.T, 0.0)

    def _kernel(self, a: NDArray[np.float64], b: NDArray[np.float64]) -> NDArray[np.float64]:
        d2 = self._pairwise_sq_dists(a, b)
        return np.exp(-d2 / (2.0 * self.epsilon**2))

    def fit(self, x: NDArray[np.float64], y: NDArray[np.float64]) -> "RBFSurrogate":
        x = np.atleast_2d(np.asarray(x, dtype=float))
        y = np.asarray(y, dtype=float).ravel()
        if x.shape[0] != y.shape[0] or x.shape[0] == 0:
            raise ValueError("x and y must have matching, non-zero length")
        xn = self._normalize(x)
        self._y_mean = float(np.mean(y))
        k = self._kernel(xn, xn)
        k += self.ridge * np.eye(k.shape[0])
        self._weights = np.linalg.solve(k, y - self._y_mean)
        self._x_train = xn
        return self

    @property
    def is_fitted(self) -> bool:
        return self._x_train is not None and self._weights is not None

    def predict(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        if not self.is_fitted:
            raise RuntimeError("surrogate is not fitted")
        xn = self._normalize(np.atleast_2d(np.asarray(x, dtype=float)))
        k = self._kernel(xn, self._x_train)  # type: ignore[arg-type]
        return self._y_mean + k @ self._weights  # type: ignore[operator]

    def uncertainty(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """Normalized distance to the nearest training point (higher = less known)."""
        if not self.is_fitted:
            raise RuntimeError("surrogate is not fitted")
        xn = self._normalize(np.atleast_2d(np.asarray(x, dtype=float)))
        d2 = self._pairwise_sq_dists(xn, self._x_train)  # type: ignore[arg-type]
        return np.sqrt(np.min(d2, axis=1))


def screen_candidates(
    surrogate: RBFSurrogate,
    candidates: NDArray[np.float64],
    *,
    keep_fraction: float = 0.5,
    explore_quantile: float = 0.75,
    min_eval: int = 1,
) -> tuple[list[int], NDArray[np.float64]]:
    """Decide which candidates to evaluate on the GPU.

    Returns ``(eval_indices, predicted_scores)``.  A candidate is evaluated if
    it is in the top ``keep_fraction`` by predicted score *or* its surrogate
    uncertainty exceeds the ``explore_quantile`` of the batch (exploration).
    At least ``min_eval`` candidates are always evaluated.
    """
    candidates = np.atleast_2d(np.asarray(candidates, dtype=float))
    n = candidates.shape[0]
    predicted = surrogate.predict(candidates)
    uncert = surrogate.uncertainty(candidates)

    n_keep = max(min_eval, int(round(keep_fraction * n)))
    order = np.argsort(-predicted)  # descending predicted score
    keep = set(int(i) for i in order[:n_keep])

    if 0.0 < explore_quantile < 1.0 and n >= 2:
        thresh = float(np.quantile(uncert, explore_quantile))
        for i in range(n):
            if uncert[i] >= thresh:
                keep.add(i)

    return sorted(keep), predicted
