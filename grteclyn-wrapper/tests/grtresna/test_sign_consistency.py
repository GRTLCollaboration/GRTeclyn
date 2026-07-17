"""Fail-closed tests for ID ↔ evolution gravitational-sign mismatch.

Catches the eval-118 / qball_traj bug: GRTresna solves with per-lump phantom
signs while grtresna_complex_scalar evolves one globally canonical field.
"""

from __future__ import annotations

import pytest

from grteclyn_wrapper.grtresna.matter.models import (
    GRTRESNA_BICOMPLEX_SCALAR_MODEL,
    GRTRESNA_COMPLEX_SCALAR_MODEL,
    GRTRESNA_INDEPENDENT_MATTER_MODEL,
)
from grteclyn_wrapper.grtresna.matter.sign_consistency import (
    SignMismatchError,
    assert_id_evolution_sign_consistency,
    check_id_evolution_sign_consistency,
    evolution_effective_signs,
    id_lump_signs,
)
from grteclyn_wrapper.grtresna.matter.wiring import (
    GRTresnaMatterMetadata,
    evolution_overrides_from_config,
)
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig


def _lump(*, exotic: int = 0, amp: float = 0.15, center=(0.0, 0.0, 0.0)) -> dict:
    return {
        "amp": amp,
        "width": 1.67,
        "center": center,
        "velocity": (0.0, 0.0, 0.0),
        "omega": 0.8,
        "mode": 0,
        "exotic": exotic,
    }


def test_id_lump_signs_from_exotic_flags() -> None:
    signs = id_lump_signs(
        [
            _lump(exotic=0),
            _lump(exotic=1, center=(1.0, 0.0, 0.0)),
            _lump(exotic=1, center=(2.0, 0.0, 0.0)),
        ]
    )
    assert signs == (1, -1, -1)


def test_evolution_effective_signs_single_complex_flattens_to_global() -> None:
    # This is the silent drop: five ID signs become five copies of +1.
    evo = evolution_effective_signs(
        GRTRESNA_COMPLEX_SCALAR_MODEL,
        scalar_sign=1,
        n_lumps=5,
    )
    assert evo == (1, 1, 1, 1, 1)


def test_catch_eval118_style_complex_scalar_exotic_mismatch() -> None:
    """Eval-118 campaign shape: mixed exotic lumps + single-complex + sign=+1."""
    cfg = GRTresnaConfig(
        matter_model=GRTRESNA_COMPLEX_SCALAR_MODEL,
        scalar_sign=1,
        scalar_mass=1.0,
        scalar_lambda=640.0,
        scalar_mu=85333.0,
        bs_omega=0.8,
        lumps=[
            _lump(exotic=0, center=(-6.0, 1.0, 2.0)),
            _lump(exotic=1, center=(3.0, 7.0, 1.0)),
            _lump(exotic=1, center=(-0.4, 1.0, -1.6)),
            _lump(exotic=1, center=(0.8, 3.0, 2.0)),
            _lump(exotic=1, center=(4.5, -0.5, -1.7)),
        ],
    )
    report = check_id_evolution_sign_consistency(cfg)
    assert report.ok is False
    assert report.id_lump_signs == (1, -1, -1, -1, -1)
    assert report.evolution_signs == (1, 1, 1, 1, 1)
    assert "bicomplex" in (report.reason or "").lower()
    with pytest.raises(SignMismatchError, match="phantom lumps"):
        assert_id_evolution_sign_consistency(cfg)


def test_catch_all_exotic_lumps_with_canonical_global_sign() -> None:
    cfg = GRTresnaConfig(
        matter_model=GRTRESNA_COMPLEX_SCALAR_MODEL,
        scalar_sign=1,
        lumps=[_lump(exotic=1), _lump(exotic=1, center=(2.0, 0.0, 0.0))],
    )
    report = check_id_evolution_sign_consistency(cfg)
    assert report.ok is False
    assert report.id_lump_signs == (-1, -1)
    assert report.evolution_signs == (1, 1)


def test_catch_mixed_signs_even_if_global_phantom() -> None:
    """Global scalar_sign=-1 still cannot represent mixed per-lump signs."""
    cfg = GRTresnaConfig(
        matter_model=GRTRESNA_COMPLEX_SCALAR_MODEL,
        scalar_sign=-1,
        lumps=[
            _lump(exotic=0),
            _lump(exotic=1, center=(2.0, 0.0, 0.0)),
        ],
    )
    report = check_id_evolution_sign_consistency(cfg)
    assert report.ok is False
    assert "mixed" in (report.reason or "").lower()


def test_ok_all_canonical_complex_scalar() -> None:
    cfg = GRTresnaConfig(
        matter_model=GRTRESNA_COMPLEX_SCALAR_MODEL,
        scalar_sign=1,
        lumps=[_lump(exotic=0), _lump(exotic=0, center=(2.0, 0.0, 0.0))],
    )
    report = assert_id_evolution_sign_consistency(cfg)
    assert report.ok is True
    assert report.id_lump_signs == (1, 1)
    assert report.evolution_signs == (1, 1)


def test_ok_all_phantom_complex_scalar_with_matching_global_sign() -> None:
    cfg = GRTresnaConfig(
        matter_model=GRTRESNA_COMPLEX_SCALAR_MODEL,
        scalar_sign=-1,
        lumps=[_lump(exotic=1), _lump(exotic=1, center=(2.0, 0.0, 0.0))],
    )
    report = assert_id_evolution_sign_consistency(cfg)
    assert report.ok is True
    assert report.evolution_signs == (-1, -1)


def test_ok_bicomplex_retains_per_lump_signs() -> None:
    cfg = GRTresnaConfig(
        matter_model=GRTRESNA_BICOMPLEX_SCALAR_MODEL,
        scalar_sign=1,
        lumps=[
            _lump(exotic=0),
            _lump(exotic=1, center=(1.0, 0.0, 0.0)),
            _lump(exotic=1, center=(2.0, 0.0, 0.0)),
        ],
    )
    report = assert_id_evolution_sign_consistency(cfg)
    assert report.ok is True
    assert report.id_lump_signs == (1, -1, -1)
    overrides = evolution_overrides_from_config(cfg)
    assert overrides["recipe_matter_model"] == GRTRESNA_BICOMPLEX_SCALAR_MODEL
    assert overrides["recipe_scalar_field_signs"] == "1 -1 -1"


def test_ok_independent_scalars_with_mixed_signs() -> None:
    cfg = GRTresnaConfig(
        matter_model=GRTRESNA_INDEPENDENT_MATTER_MODEL,
        lumps=[
            _lump(exotic=0),
            _lump(exotic=1, center=(2.0, 0.0, 0.0)),
        ],
    )
    report = assert_id_evolution_sign_consistency(cfg)
    assert report.ok is True
    overrides = evolution_overrides_from_config(cfg)
    assert overrides["recipe_scalar_field_signs"] == "1 -1"


def test_metadata_from_complex_scalar_drops_per_lump_signs() -> None:
    """Freeze the wiring bug the consistency rail must catch."""
    cfg = GRTresnaConfig(
        matter_model=GRTRESNA_COMPLEX_SCALAR_MODEL,
        scalar_sign=1,
        lumps=[_lump(exotic=0), _lump(exotic=1, center=(2.0, 0.0, 0.0))],
    )
    meta = GRTresnaMatterMetadata.from_config(cfg)
    assert meta.scalar_field_signs == ()
    assert meta.scalar_sign == 1
    # Metadata alone looks "canonical"; the rail must inspect the lumps.
    assert check_id_evolution_sign_consistency(cfg).ok is False
