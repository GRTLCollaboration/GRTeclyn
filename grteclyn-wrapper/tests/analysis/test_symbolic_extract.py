"""Tests for symbolic metric extraction and cross-code verification hooks."""

from __future__ import annotations

from grteclyn_wrapper.analysis.symbolic_extract import (
    SymbolicMetricForm,
    verify_analytic_form,
)


def test_verify_analytic_form_accepts_good_fit():
    form = SymbolicMetricForm(
        expressions={"chi": "1.0", "lapse": "1.0", "shift1": "0.0"},
        slice="xz_midplane_axisymmetric",
        r2_scores={"chi": 0.95, "lapse": 0.9, "shift1": 0.85},
    )
    payload = verify_analytic_form(form, geodesic_f_geo=0.05)
    assert payload["analytic_form"] is True
    assert payload["resolution_converged"] is True
