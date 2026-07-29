"""The post-load gate must read the whole-hierarchy constraint columns.

Columns 2-3 of ``constraint_norms.dat`` are computed on level 0 only, over
every level-0 cell including the ones hidden under the refinement, and averaged
across a domain that is overwhelmingly empty.  They are the right input for the
pump governor and the wrong number to gate initial data on.  Columns 12/16/17
are the composite: all levels, covered cells dropped, plus an undiluted peak
and a refined-region-only norm.

The trap these tests exist to hold shut: on a single-level run column 17 is
written as exactly ``0.0``.  That is *unmeasured*, and it looks like perfect.
"""

from __future__ import annotations

from pathlib import Path

from grteclyn_wrapper.metrics.diagnostics.constraints import read_constraint_metrics
from grteclyn_wrapper.projection.postload_gate import (
    PostLoadGateConfig,
    evaluate_constraint_gate,
)

# time, L2_Ham, L2_Mom, min_rho, max_rho, int_neg_rho, L2_Ham_rel, L2_Mom_rel,
# pump_force_L2, governor, pump_fi_L2, L2_Ham_amr, L2_Mom_amr, L2_Ham_amr_rel,
# L2_Mom_amr_rel, Linf_Ham_amr, L2_Ham_amr_ref
_HEADER = (
    "# time L2_Ham L2_Mom min_rho_req max_rho_req integral_neg_rho L2_Ham_rel "
    "L2_Mom_rel pump_force_L2 governor pump_fi_L2 L2_Ham_amr L2_Mom_amr "
    "L2_Ham_amr_rel L2_Mom_amr_rel Linf_Ham_amr L2_Ham_amr_ref\n"
)


def _write(path: Path, *rows: list[float], header: bool = True) -> Path:
    text = _HEADER if header else ""
    for row in rows:
        text += " ".join(f"{v:.10e}" for v in row) + "\n"
    path.write_text(text, encoding="utf-8")
    return path


def _wide_row(
    *,
    l2_ham: float,
    l2_ham_amr: float,
    linf_ham_amr: float,
    l2_ham_amr_ref: float,
) -> list[float]:
    row = [0.01, l2_ham, 1.0e-4, -1.0e-3, 1.0e-3, 0.0, 0.1, 0.1, 0.0, 1.0, 0.0]
    row += [l2_ham_amr, 1.0e-4, 0.1, 0.1, linf_ham_amr, l2_ham_amr_ref]
    return row


def test_the_composite_columns_are_read(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "constraint_norms.dat",
        _wide_row(
            l2_ham=1.0e-3,
            l2_ham_amr=4.0e-3,
            linf_ham_amr=2.5e-1,
            l2_ham_amr_ref=7.0e-3,
        ),
    )
    metrics = read_constraint_metrics(path)
    assert metrics is not None
    assert metrics.max_hamiltonian_l2_amr == 4.0e-3
    assert metrics.max_hamiltonian_linf_amr == 2.5e-1
    assert metrics.max_hamiltonian_l2_amr_ref == 7.0e-3
    assert metrics.has_refined_region is True


def test_a_legacy_file_reports_no_composite_columns(tmp_path: Path) -> None:
    """Campaigns A-F wrote 11 columns; those runs must still parse."""
    path = tmp_path / "constraint_norms.dat"
    path.write_text(
        "# time L2_Ham L2_Mom min_rho_req max_rho_req integral_neg_rho\n"
        "1.0000000e-02 1.3458257985e-02 3.3270691728e-04 "
        "-3.3605017576e-03 2.9936897853e-03 1.9440106788e+01\n",
        encoding="utf-8",
    )
    metrics = read_constraint_metrics(path)
    assert metrics is not None
    assert metrics.max_hamiltonian_l2 == 1.3458257985e-02
    assert metrics.max_hamiltonian_l2_amr is None
    assert metrics.has_refined_region is False

    gate = evaluate_constraint_gate(path, config=PostLoadGateConfig())
    assert any("predates" in note for note in gate.notes)


def test_a_single_level_run_does_not_pass_on_an_empty_refined_region(
    tmp_path: Path,
) -> None:
    """Col 17 is 0.0 when no level 1 exists.  That must not read as clean."""
    path = _write(
        tmp_path / "constraint_norms.dat",
        _wide_row(
            l2_ham=1.0e-3,
            l2_ham_amr=4.0e-3,
            linf_ham_amr=2.5e-1,
            l2_ham_amr_ref=0.0,
        ),
    )
    metrics = read_constraint_metrics(path)
    assert metrics is not None
    assert metrics.max_hamiltonian_l2_amr_ref == 0.0
    assert metrics.has_refined_region is False

    gate = evaluate_constraint_gate(path, config=PostLoadGateConfig())
    assert any("unmeasured" in note for note in gate.notes)
    assert not any("refined-region L2" in note for note in gate.notes)


def test_the_composite_peak_can_fail_data_the_domain_mean_accepts(
    tmp_path: Path,
) -> None:
    """The whole point: a violation localised on the lump is invisible to col 2.

    ``L2_Ham`` is a hair under the legacy 1e-2 threshold because the empty
    domain dilutes it, while the undiluted peak is two orders of magnitude
    larger.  Before this change the gate saw only the former and passed.
    """
    path = _write(
        tmp_path / "constraint_norms.dat",
        _wide_row(
            l2_ham=9.0e-3,
            l2_ham_amr=8.0e-3,
            linf_ham_amr=3.0,
            l2_ham_amr_ref=5.0e-2,
        ),
    )

    reporting_only = evaluate_constraint_gate(path, config=PostLoadGateConfig())
    assert reporting_only.passed is True
    assert reporting_only.max_hamiltonian_linf_amr == 3.0
    assert any("no composite threshold configured" in n for n in reporting_only.notes)

    gated = evaluate_constraint_gate(
        path,
        config=PostLoadGateConfig(max_hamiltonian_linf_amr=1.0e-1),
    )
    assert gated.passed is False
    assert gated.reason is not None
    assert "composite" in gated.reason


def test_the_composite_values_reach_the_result(tmp_path: Path) -> None:
    """A threshold can only be calibrated from numbers that get recorded."""
    path = _write(
        tmp_path / "constraint_norms.dat",
        _wide_row(
            l2_ham=1.0e-3,
            l2_ham_amr=4.0e-3,
            linf_ham_amr=2.5e-1,
            l2_ham_amr_ref=7.0e-3,
        ),
    )
    gate = evaluate_constraint_gate(path, config=PostLoadGateConfig())
    assert gate.max_hamiltonian_l2_amr == 4.0e-3
    assert gate.max_hamiltonian_linf_amr == 2.5e-1
    assert gate.max_hamiltonian_l2_amr_ref == 7.0e-3
    assert gate.has_refined_region is True
