"""Shared helpers for diagnostic metric parsers."""

from __future__ import annotations


def positive_fractional_change(final: float, initial: float, *, floor: float) -> float:
    return max(0.0, final - initial) / max(abs(initial), floor)


def positive_fractional_drop(initial: float, minimum: float, *, floor: float) -> float:
    return max(0.0, initial - minimum) / max(abs(initial), floor)
