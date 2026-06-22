"""Tests for central plotfile extraction helpers."""

from __future__ import annotations

from grteclyn_wrapper.metrics.diagnostics.central import CENTRAL_TIMESERIES_HEADER
from grteclyn_wrapper.visualisation.process_wave.consume_plotfiles.extraction.central import (
    _extract_central_timeseries_line,
    _field_at_point,
    _read_scalar,
    _scalar_activity,
)


class _FakePoint:
    def __init__(self, values: dict[tuple[str, str], float]):
        self._values = values

    def __getitem__(self, key):
        return self._values[key]


class _FakeDataset:
    def __init__(self, fields: dict[str, float]):
        self.field_list = [("boxlib", name) for name in fields]
        self._fields = fields

    def point(self, _center):
        return _FakePoint({("boxlib", k): v for k, v in self._fields.items()})


def test_central_timeseries_header_format() -> None:
    assert "rho_req" in CENTRAL_TIMESERIES_HEADER
    assert "scalar_activity" in CENTRAL_TIMESERIES_HEADER
    # Geometric (spacetime-splash) columns must be present.
    assert "chi" in CENTRAL_TIMESERIES_HEADER
    assert "K" in CENTRAL_TIMESERIES_HEADER
    assert "weyl4" in CENTRAL_TIMESERIES_HEADER


def test_field_at_point_resolves_weyl4_re() -> None:
    """weyl4 resolves to the Weyl4_Re plotfile field."""
    ds = _FakeDataset({"Weyl4_Re": -0.013})
    point = ds.point([0.0, 0.0, 0.0])
    assert _field_at_point(ds, point, "weyl4") == -0.013


def test_read_scalar_finite() -> None:
    assert _read_scalar([1.25]) == 1.25


def test_field_at_point_returns_value() -> None:
    ds = _FakeDataset({"rho_req": 0.42})
    point = ds.point([0.0, 0.0, 0.0])
    assert _field_at_point(ds, point, "rho_req") == 0.42


def test_scalar_activity_uses_direct_field() -> None:
    ds = _FakeDataset({"scalar_activity": 0.17})
    point = ds.point([0.0, 0.0, 0.0])
    assert _scalar_activity(ds, point) == 0.17


def test_scalar_activity_fallback_from_components() -> None:
    ds = _FakeDataset({"phi": 3.0, "Pi": 4.0})
    point = ds.point([0.0, 0.0, 0.0])
    assert _scalar_activity(ds, point) == 5.0


def test_extract_line_format(monkeypatch) -> None:
    fields = {
        "rho_req": 1.0e-3,
        "lapse": 0.95,
        "scalar_activity": 0.02,
        "phi": 0.1,
        "phi_lump0": 0.0,
    }

    class _FakeYt:
        @staticmethod
        def load(_path):
            return _FakeDataset(fields)

    monkeypatch.setitem(
        __import__("sys").modules,
        "yt",
        _FakeYt(),
    )
    line = _extract_central_timeseries_line(
        "/fake/plt0000",
        t=1.5,
        center=(0.0, 0.0, 0.0),
    )
    assert line is not None
    assert line.startswith("1.5")
    assert "1.0000000000000000e-03" in line


def test_extract_line_includes_geometric_fields(monkeypatch) -> None:
    """chi, K, Weyl4_Re are extracted into the timeseries line (13 columns)."""
    fields = {
        "rho_req": 1.0e-3,
        "lapse": 0.95,
        "scalar_activity": 0.02,
        "phi": 0.1,
        "phi_lump0": 0.0,
        "chi": 0.42,
        "K": 1.7,
        "Weyl4_Re": -2.5e-3,
    }

    class _FakeYt:
        @staticmethod
        def load(_path):
            return _FakeDataset(fields)

    monkeypatch.setitem(__import__("sys").modules, "yt", _FakeYt())
    line = _extract_central_timeseries_line(
        "/fake/plt0000",
        t=1.5,
        center=(0.0, 0.0, 0.0),
    )
    assert line is not None
    cols = line.split()
    # 13 columns: time + 9 original + chi + K + weyl4
    assert len(cols) == 13
    assert float(cols[10]) == 0.42      # chi
    assert float(cols[11]) == 1.7       # K
    assert float(cols[12]) == -2.5e-3   # weyl4
