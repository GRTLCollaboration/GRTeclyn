"""Per-sector barycentre stream (the Bondi-dipole runaway diagnostic).

The one property the aggregate confinement stream cannot provide: for a
mixed-sign pair the total-matter barycentre cancels while each sector's own
centroid drifts.  These tests pin the split-and-centroid behaviour, including
the two cases the field-sniff fallback gets wrong (all-phantom / all-canonical
runs), which is why the extractor takes the model tag explicitly.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

yt = pytest.importorskip("yt")

from grteclyn_wrapper.visualisation.process_wave.consume_plotfiles.extraction.sector_barycenters import (
    SECTOR_BARYCENTERS_COLUMNS,
    _extract_sector_barycenters_line,
)

_BICOMPLEX = "grtresna_bicomplex_scalar"


def _parse(row: str | None) -> dict[str, float]:
    assert row is not None
    cols = [float(x) for x in row.split()]
    assert len(cols) == len(SECTOR_BARYCENTERS_COLUMNS)
    return dict(zip(SECTOR_BARYCENTERS_COLUMNS, cols))


def _pair_ds(*, canon_center: float = 4.0, phantom_center: float = 12.0):
    """Uniform grid holding one canonical and one phantom blob on the x axis.

    Bicomplex layout (StateVariables.hpp): phi/Pi = canonical Re/momentum,
    phi_lump0/Pi_lump0 = canonical Im, phi_lump1..2 = phantom Re/Im.
    """
    n = 16
    zeros = np.zeros((n, n, n), dtype=np.float64)
    xs = np.arange(n) + 0.5  # cell centres, dx = 1

    def blob(center_x: float) -> np.ndarray:
        w = np.zeros((n, n, n), dtype=np.float64)
        w[np.abs(xs - center_x) < 1.0, :, :] = 1.0
        return w

    data = {
        ("stream", "phi"): blob(canon_center),
        ("stream", "Pi"): zeros.copy(),
        ("stream", "phi_lump0"): zeros.copy(),
        ("stream", "Pi_lump0"): zeros.copy(),
        ("stream", "phi_lump1"): blob(phantom_center),
        ("stream", "Pi_lump1"): zeros.copy(),
        ("stream", "phi_lump2"): zeros.copy(),
        ("stream", "Pi_lump2"): zeros.copy(),
    }
    return yt.load_uniform_grid(
        data, (n, n, n),
        bbox=np.array([[0.0, 16.0], [0.0, 16.0], [0.0, 16.0]]),
    )


def _extract(ds, **kwargs) -> str | None:
    orig = yt.load
    yt.load = lambda *_a, **_k: ds
    try:
        return _extract_sector_barycenters_line("x", t=1.5, **kwargs)
    finally:
        yt.load = orig


def test_pair_run_reports_each_sector_centroid():
    row = _parse(_extract(_pair_ds(), matter_model=_BICOMPLEX))
    assert row["time"] == pytest.approx(1.5)
    assert row["bary_x_canon"] == pytest.approx(4.0, abs=0.01)
    assert row["bary_x_phantom"] == pytest.approx(12.0, abs=0.01)
    # Blobs are slabs: centred in y/z at the domain midpoint.
    assert row["bary_y_canon"] == pytest.approx(8.0, abs=0.01)
    assert row["bary_z_phantom"] == pytest.approx(8.0, abs=0.01)
    assert row["total_canon"] > 0.0
    assert row["total_phantom"] > 0.0


def test_sector_drift_is_visible_where_aggregate_cancels():
    """Move both blobs +2 in x: each sector centroid moves, symmetrically."""
    before = _parse(_extract(_pair_ds(canon_center=4.0, phantom_center=12.0), matter_model=_BICOMPLEX))
    after = _parse(_extract(_pair_ds(canon_center=6.0, phantom_center=14.0), matter_model=_BICOMPLEX))
    assert after["bary_x_canon"] - before["bary_x_canon"] == pytest.approx(2.0, abs=0.02)
    assert after["bary_x_phantom"] - before["bary_x_phantom"] == pytest.approx(2.0, abs=0.02)


def test_empty_sector_reports_nan_not_origin():
    """All-phantom run: canonical centroid must be nan, never (0,0,0).

    This is exactly the configuration the field sniff cannot classify (zero
    canonical looks like independent scalars) -- the explicit tag decides.
    """
    n = 8
    zeros = np.zeros((n, n, n), dtype=np.float64)
    blob = zeros.copy()
    blob[2:4, 2:4, 2:4] = 1.0
    data = {
        ("stream", "phi"): zeros.copy(),
        ("stream", "Pi"): zeros.copy(),
        ("stream", "phi_lump0"): zeros.copy(),
        ("stream", "Pi_lump0"): zeros.copy(),
        ("stream", "phi_lump1"): blob,
        ("stream", "Pi_lump1"): zeros.copy(),
        ("stream", "phi_lump2"): zeros.copy(),
        ("stream", "Pi_lump2"): zeros.copy(),
    }
    ds = yt.load_uniform_grid(
        data, (n, n, n),
        bbox=np.array([[0.0, 8.0], [0.0, 8.0], [0.0, 8.0]]),
    )
    row = _parse(_extract(ds, matter_model=_BICOMPLEX))
    assert math.isnan(row["bary_x_canon"])
    assert math.isnan(row["rms_radius_canon"])
    assert row["total_phantom"] > 0.0
    assert not math.isnan(row["bary_x_phantom"])


def test_single_sector_model_yields_no_row():
    """A model with no canonical/phantom split writes nothing, not nan spam."""
    n = 8
    data = {
        ("stream", "phi"): np.full((n, n, n), 0.1, dtype=np.float64),
        ("stream", "Pi"): np.zeros((n, n, n), dtype=np.float64),
    }
    ds = yt.load_uniform_grid(
        data, (n, n, n),
        bbox=np.array([[0.0, 8.0], [0.0, 8.0], [0.0, 8.0]]),
    )
    assert _extract(ds, matter_model="") is None
