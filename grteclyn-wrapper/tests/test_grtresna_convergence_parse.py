import math
from pathlib import Path

from grteclyn_wrapper.grtresna.solver import parse_convergence


def _write_errors(tmp_path: Path, rows: list[str]) -> Path:
    work = tmp_path / "grtresna"
    work.mkdir()
    (work / "Ham_and_Mom_errors.txt").write_text(
        "iteration ham_pct mom_pct\n" + "\n".join(rows) + "\n"
    )
    return work


def test_static_matter_nan_momentum_is_treated_as_satisfied(tmp_path):
    # Static (momentum-free) matter: relative momentum error is 0/0 = -nan even
    # though the Hamiltonian converged.  Must be read as a satisfied (0%)
    # momentum constraint, not a nonfinite rejection.
    work = _write_errors(tmp_path, ["48 7.59e-03 -nan", "49 7.57e-03 -nan"])
    conv = parse_convergence(work)
    assert conv is not None
    assert math.isfinite(conv["ham_pct"]) and conv["ham_pct"] < 0.05
    assert conv["mom_pct"] == 0.0


def test_full_divergence_stays_nonfinite(tmp_path):
    # A genuine solver divergence shows up as NaN in *both* Ham and Mom and must
    # remain nonfinite so the convergence gate still rejects it.
    work = _write_errors(tmp_path, ["49 -nan -nan"])
    conv = parse_convergence(work)
    assert conv is not None
    assert not math.isfinite(conv["ham_pct"])
    assert not math.isfinite(conv["mom_pct"])


def test_moving_matter_finite_residuals_untouched(tmp_path):
    work = _write_errors(tmp_path, ["49 0.012 0.031"])
    conv = parse_convergence(work)
    assert conv == {"iteration": 49, "ham_pct": 0.012, "mom_pct": 0.031}


def test_ham_divergence_with_nan_mom_still_flags_ham(tmp_path):
    # Mom coerced to 0, but a diverged Hamiltonian must still be visible so the
    # gate rejects on Ham.
    work = _write_errors(tmp_path, ["49 100.0 -nan"])
    conv = parse_convergence(work)
    assert conv["ham_pct"] == 100.0
    assert conv["mom_pct"] == 0.0
