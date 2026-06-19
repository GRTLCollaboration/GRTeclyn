"""Unit tests for gridinit export window / spacing policy."""

from __future__ import annotations

import numpy as np
import pytest

from grteclyn_wrapper.grtresna import io
from grteclyn_wrapper.grtresna.domain import GRTresnaDomainConfig


def test_target_span_slice_respects_source_origin() -> None:
    dx_target = np.array([0.5, 0.5, 0.5])
    # Chombo cell at x in [64, 66) with dx_lev=2 should map near evolution centre
    # when the export window starts at source_origin=32.
    span = io._target_span_slice(
        32, 0, 2.0, dx_target, 128, target_origin=32.0,
    )
    assert span.start >= 62
    assert span.stop <= 70


def test_half_z_domain_keeps_legacy_full_solve_export() -> None:
    domain = GRTresnaDomainConfig(full_z=False, l_full=80.0, n_full=96, grtresna_nx=48)
    export = domain.gridinit_export_spec()
    assert export.lx == export.ly == export.lz == 128.0
    assert export.source_origin == (0.0, 0.0, 0.0)
