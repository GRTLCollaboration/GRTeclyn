"""Serialization helpers for frozen metric dataclasses."""

from __future__ import annotations

import numpy as np


def dataclass_to_dict(value: object) -> object:
    if hasattr(value, "__dataclass_fields__"):
        return {
            key: dataclass_to_dict(getattr(value, key))
            for key in value.__dataclass_fields__
        }
    if isinstance(value, dict):
        return {key: dataclass_to_dict(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [dataclass_to_dict(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value
