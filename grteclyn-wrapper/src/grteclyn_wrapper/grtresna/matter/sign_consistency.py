"""Fail-closed ID ↔ evolution gravitational-sign consistency rail.

GRTresna's boson-star constraint solve can apply *per-lump* exotic signs to
``T_ab``.  The single-field ``grtresna_complex_scalar`` evolution path only
supports one global ``recipe_scalar_sign``.  Wiring exotic lumps into that
path silently drops the per-lump signs (eval-118 / qball_traj bug class).

This module detects that mismatch before GPU time is spent.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from ..fields.lump import lump_sign
from .models import (
    is_bicomplex_scalar_model,
    is_complex_scalar_model,
)

logger = logging.getLogger(__name__)

# Explicit opt-out for legacy debugging only. Default is fail-closed.
_ALLOW_MISMATCH_ENV = "GRTRESNA_ALLOW_SIGN_MISMATCH"


@dataclass(frozen=True)
class SignConsistencyReport:
    """Result of an ID/evolution sign audit."""

    ok: bool
    matter_model: str
    id_lump_signs: tuple[int, ...]
    evolution_signs: tuple[int, ...]
    reason: str | None = None

    def raise_if_mismatch(self) -> None:
        if not self.ok:
            raise SignMismatchError(self.reason or "ID/evolution sign mismatch")


class SignMismatchError(ValueError):
    """Raised when constraint-solve signs disagree with evolution ``T_ab``."""


def _active_lumps(lumps: Sequence[Mapping[str, Any]] | None) -> list[Mapping[str, Any]]:
    out: list[Mapping[str, Any]] = []
    for lump in lumps or ():
        try:
            if float(lump.get("amp", 0.0)) == 0.0:
                continue
        except (TypeError, ValueError):
            continue
        out.append(lump)
    return out


def id_lump_signs(lumps: Sequence[Mapping[str, Any]] | None) -> tuple[int, ...]:
    """Per-lump gravitational signs used by the GRTresna constraint solve."""
    return tuple(lump_sign(lump) for lump in _active_lumps(lumps))


def evolution_effective_signs(
    matter_model: str,
    *,
    scalar_sign: int = 1,
    scalar_field_signs: Sequence[int] | None = None,
    n_lumps: int = 0,
) -> tuple[int, ...]:
    """Signs that GRTeclyn will apply to ``T_ab`` for this matter model."""
    if is_bicomplex_scalar_model(matter_model):
        return tuple(int(s) for s in (scalar_field_signs or ()))
    if is_complex_scalar_model(matter_model):
        # One summed field: every lump is forced to the global sign.
        sign = 1 if int(scalar_sign) >= 0 else -1
        return tuple(sign for _ in range(max(0, int(n_lumps))))
    # Independent scalars: per-channel signs travel through recipe_scalar_field_signs.
    return tuple(int(s) for s in (scalar_field_signs or ()))


def check_id_evolution_sign_consistency(
    cfg: Any,
) -> SignConsistencyReport:
    """Return ok=False when ID exotic signs are not retained in evolution ``T_ab``.

    Catches the eval-118 failure mode:
    ``grtresna_complex_scalar`` + per-lump ``exotic=1`` + global ``scalar_sign=+1``.
    """
    matter_model = str(getattr(cfg, "matter_model", "") or "")
    lumps = list(getattr(cfg, "lumps", ()) or ())
    id_signs = id_lump_signs(lumps)
    scalar_sign = int(getattr(cfg, "scalar_sign", 1) or 1)

    if is_bicomplex_scalar_model(matter_model):
        evo_signs = evolution_effective_signs(
            matter_model,
            scalar_field_signs=id_signs,
            n_lumps=len(id_signs),
        )
        if id_signs != evo_signs:
            return SignConsistencyReport(
                ok=False,
                matter_model=matter_model,
                id_lump_signs=id_signs,
                evolution_signs=evo_signs,
                reason=(
                    "bicomplex field signs do not match ID lump signs: "
                    f"id={id_signs} evo={evo_signs}"
                ),
            )
        return SignConsistencyReport(
            ok=True,
            matter_model=matter_model,
            id_lump_signs=id_signs,
            evolution_signs=evo_signs,
        )

    if is_complex_scalar_model(matter_model):
        evo_signs = evolution_effective_signs(
            matter_model,
            scalar_sign=scalar_sign,
            n_lumps=len(id_signs),
        )
        if not id_signs:
            return SignConsistencyReport(
                ok=True,
                matter_model=matter_model,
                id_lump_signs=id_signs,
                evolution_signs=evo_signs,
            )
        if any(s < 0 for s in id_signs) and scalar_sign >= 0:
            return SignConsistencyReport(
                ok=False,
                matter_model=matter_model,
                id_lump_signs=id_signs,
                evolution_signs=evo_signs,
                reason=(
                    "single-complex evolution uses global canonical sign +1, but "
                    f"GRTresna ID has phantom lumps {id_signs}; use "
                    "grtresna_bicomplex_scalar to retain per-lump signs"
                ),
            )
        if len(set(id_signs)) > 1:
            return SignConsistencyReport(
                ok=False,
                matter_model=matter_model,
                id_lump_signs=id_signs,
                evolution_signs=evo_signs,
                reason=(
                    "mixed canonical/phantom lumps cannot be represented by "
                    f"grtresna_complex_scalar (id={id_signs}); use "
                    "grtresna_bicomplex_scalar"
                ),
            )
        if id_signs != evo_signs:
            return SignConsistencyReport(
                ok=False,
                matter_model=matter_model,
                id_lump_signs=id_signs,
                evolution_signs=evo_signs,
                reason=(
                    "ID lump signs disagree with single-complex evolution sign: "
                    f"id={id_signs} evo={evo_signs}"
                ),
            )
        return SignConsistencyReport(
            ok=True,
            matter_model=matter_model,
            id_lump_signs=id_signs,
            evolution_signs=evo_signs,
        )

    # Independent scalars: evolution carries recipe_scalar_field_signs.
    evo_signs = evolution_effective_signs(
        matter_model,
        scalar_field_signs=id_signs,
        n_lumps=len(id_signs),
    )
    ok = id_signs == evo_signs
    return SignConsistencyReport(
        ok=ok,
        matter_model=matter_model,
        id_lump_signs=id_signs,
        evolution_signs=evo_signs,
        reason=None
        if ok
        else (
            "independent-scalar evolution signs disagree with ID: "
            f"id={id_signs} evo={evo_signs}"
        ),
    )


def allow_sign_mismatch() -> bool:
    """True only when ``GRTRESNA_ALLOW_SIGN_MISMATCH=1`` (legacy debug)."""
    return os.environ.get(_ALLOW_MISMATCH_ENV, "0").strip().lower() in (
        "1",
        "true",
        "yes",
    )


def assert_id_evolution_sign_consistency(cfg: Any) -> SignConsistencyReport:
    """Raise :class:`SignMismatchError` when ID/evolution signs disagree.

    Honours ``GRTRESNA_ALLOW_SIGN_MISMATCH=1`` as a warning-only downgrade for
    debugging corrupted archives. Production / QD default is hard-fail.
    """
    report = check_id_evolution_sign_consistency(cfg)
    if report.ok:
        return report
    if allow_sign_mismatch():
        logger.warning(
            "ID/evolution sign mismatch ALLOWED by %s=1: %s",
            _ALLOW_MISMATCH_ENV,
            report.reason,
        )
        return report
    report.raise_if_mismatch()
    return report


def enforce_id_evolution_sign_consistency(cfg: Any) -> SignConsistencyReport:
    """Alias used at config-build / solve entry points (fail-closed)."""
    return assert_id_evolution_sign_consistency(cfg)
