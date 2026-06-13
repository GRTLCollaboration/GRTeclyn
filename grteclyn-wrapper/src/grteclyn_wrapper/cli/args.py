"""Shared argparse value parsing helpers."""

from __future__ import annotations

import argparse
from typing import Any


def parse_override(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError(f"Override must be key=value, got {value!r}")
    key, raw = value.split("=", 1)
    key = key.strip()
    raw = raw.strip()
    if not key:
        raise argparse.ArgumentTypeError("Override key cannot be empty")
    return key, raw


def coerce_value(raw: str) -> Any:
    for caster in (int, float):
        try:
            return caster(raw)
        except ValueError:
            pass
    if raw.lower() in {"true", "false"}:
        return raw.lower() == "true"
    return raw


def collect_overrides(pairs: list[tuple[str, str]]) -> dict[str, Any]:
    return {key: coerce_value(value) for key, value in pairs}


def parse_score_weights(pairs: list[tuple[str, str]]) -> dict[str, float]:
    weights: dict[str, float] = {}
    for key, raw in pairs:
        weights[key] = float(raw)
    return weights
