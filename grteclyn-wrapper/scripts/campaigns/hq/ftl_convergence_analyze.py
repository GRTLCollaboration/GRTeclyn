#!/usr/bin/env python3
"""Generic FTL validation / Richardson analyzer for HQ RadialRecipe runs.

Unlike the rotating-wormhole hard-coded postrun script, this tool accepts
explicit run paths and grid spacings, extracts f_geo / confinement /
constraints / Psi4 summaries, and reports unequal-spacing Richardson results
only when monotonic. Non-monotonic series are marked ``rejected_extrapolation``.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

GAUGE_SPEED = math.sqrt(2.0)  # 1+log gauge modes ~ √2 c


def unequal_richardson(
    values: dict[str, float],
    spacings: dict[str, float],
    *,
    coarse: str = "coarse",
    medium: str = "medium",
    fine: str = "fine",
) -> dict[str, float | None | str]:
    """Unequal-spacing Richardson order / continuum / fine-grid GCI (%).

    Returns status ``ok`` or ``rejected_extrapolation`` (non-monotonic / no root).
    """
    required = (coarse, medium, fine)
    for key in required:
        if key not in values or key not in spacings:
            raise KeyError(f"missing key {key!r} in values/spacings")
    qc, qm, qf = (float(values[k]) for k in required)
    hc, hm, hf = (float(spacings[k]) for k in required)
    result: dict[str, float | None | str] = {
        coarse: qc,
        medium: qm,
        fine: qf,
        "order": None,
        "continuum": None,
        "gci_fine_percent": None,
        "status": "rejected_extrapolation",
        "reason": None,
    }
    dcm, dmf = qc - qm, qm - qf
    if dcm * dmf <= 0.0 or abs(dmf) < 1.0e-15:
        result["reason"] = "non_monotonic_or_zero_dmf"
        return result

    target = dcm / dmf

    def residual(order: float) -> float:
        return ((hc**order - hm**order) / (hm**order - hf**order)) - target

    lo, hi = 1.0e-6, 20.0
    flo, fhi = residual(lo), residual(hi)
    if flo * fhi > 0.0:
        result["reason"] = "no_positive_order_root"
        return result
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        fmid = residual(mid)
        if flo * fmid <= 0.0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    order = 0.5 * (lo + hi)
    ratio = (hm / hf) ** order
    continuum = qf + (qf - qm) / (ratio - 1.0)
    result["order"] = order
    result["continuum"] = continuum
    if abs(qf) > 0.0:
        result["gci_fine_percent"] = 125.0 * abs(qf - qm) / (abs(qf) * (ratio - 1.0))
    result["status"] = "ok"
    result["reason"] = None
    return result


def relative_diff(a: float, b: float) -> float | None:
    denom = max(abs(a), abs(b), 1.0e-15)
    return abs(a - b) / denom


def trusted_time(
    *,
    box_half: float,
    source_radius: float = 0.0,
    extraction_radius: float = 0.0,
    gauge_speed: float = GAUGE_SPEED,
) -> float:
    """Earliest gauge-mode contamination time from outer boundary (√2 c)."""
    # Distance from source/extraction to nearest outer face ~ L/2 - r.
    r_obs = max(source_radius, extraction_radius)
    travel = max(0.0, box_half - r_obs)
    return travel / gauge_speed


def _load_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _score_body(score_doc: dict[str, Any] | None) -> dict[str, Any]:
    if not score_doc:
        return {}
    body = score_doc.get("score")
    return body if isinstance(body, dict) else score_doc


def extract_run_summary(run_dir: Path) -> dict[str, Any]:
    """Extract validation observables from one HQ episode directory."""
    run_dir = run_dir.expanduser().resolve()
    geo = _load_json(run_dir / "small_data" / "evolving_geodesic.json") or {}
    score = _score_body(_load_json(run_dir / "score.json"))
    comps = score.get("components") or {}
    meta = _load_json(run_dir / "metadata.json") or {}
    overrides = meta.get("overrides") or {}

    constraints_path = run_dir / "data" / "constraint_norms.dat"
    ham_max = mom_max = None
    if constraints_path.is_file():
        data = np.loadtxt(constraints_path)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.size and data.shape[1] >= 3:
            ham_max = float(np.max(np.abs(data[:, 1])))
            mom_max = float(np.max(np.abs(data[:, 2])))

    psi4_path = run_dir / "small_data" / "psi4_mode_l2m0.dat"
    psi4_info: dict[str, Any] = {"present": psi4_path.is_file(), "n_rows": 0}
    if psi4_path.is_file():
        rows = [
            ln
            for ln in psi4_path.read_text(encoding="utf-8").splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")
        ]
        psi4_info["n_rows"] = len(rows)
        # Peak |r Ψ₄| at first Re/Im pair after time (R=first radius in file).
        if rows:
            arr = np.loadtxt(psi4_path)
            if arr.ndim == 1:
                arr = arr.reshape(1, -1)
            if arr.shape[1] >= 3:
                re_im = arr[:, 1:3]
                amp = np.sqrt(re_im[:, 0] ** 2 + re_im[:, 1] ** 2)
                psi4_info["peak_abs_rpsi4_first_radius"] = float(np.max(amp))

    emit_sweep = geo.get("emit_sweep") or []
    l_full = float(overrides.get("L_full") or 0.0)
    n_full = int(overrides.get("N_full") or 0)
    h = (l_full / n_full) if n_full else None

    summary = {
        "path": str(run_dir),
        "exists": run_dir.is_dir(),
        "git_commit": meta.get("git_commit"),
        "simulation_exit_code": meta.get("simulation_exit_code"),
        "L_full": l_full or None,
        "N_full": n_full or None,
        "h": h,
        "stop_time": overrides.get("stop_time"),
        "max_level": overrides.get("max_level"),
        "f_geo": geo.get("f_geo"),
        "f_geo_frozen_peak": geo.get("f_geo_frozen_peak"),
        "t_emit": geo.get("t_emit"),
        "n_rays": geo.get("n_rays"),
        "n_reached": geo.get("n_reached"),
        "h_quality_ok": geo.get("h_quality_ok"),
        "max_h_rel_drift": geo.get("max_h_rel_drift"),
        "emit_sweep": emit_sweep,
        "emit_sweep_n": len(emit_sweep),
        "score_total": score.get("total"),
        "structural_persistence": comps.get("structural_persistence"),
        "confinement_final_frac": comps.get("confinement_final_frac"),
        "confinement_spread_ratio": comps.get("confinement_spread_ratio"),
        "numerical_survival": comps.get("numerical_survival"),
        "exotic_penalty": comps.get("exotic_penalty"),
        "ftl_geo_evolving": comps.get("ftl_geo_evolving"),
        "ham_max": ham_max,
        "mom_max": mom_max,
        "psi4": psi4_info,
        "completeness": {
            "geodesic": bool(geo),
            "emit_sweep": len(emit_sweep) >= 2,
            "score": bool(score),
            "constraints": constraints_path.is_file(),
            "psi4_l2m0": psi4_path.is_file() and psi4_info["n_rows"] > 0,
            "psi4_l2_all": (run_dir / "small_data" / "psi4_mode_l2_all.dat").is_file(),
            "psi4_directional": (run_dir / "small_data" / "psi4_directional.dat").is_file(),
        },
    }
    if l_full:
        half = 0.5 * l_full
        summary["trusted_time_source"] = trusted_time(box_half=half, source_radius=0.0)
        summary["trusted_time_r24"] = trusted_time(box_half=half, extraction_radius=24.0)
    return summary


def assess_completeness_for_validation(summary: dict[str, Any]) -> list[str]:
    """Return list of missing requirements for full Phase-2+ validation."""
    missing: list[str] = []
    c = summary.get("completeness") or {}
    if not c.get("geodesic"):
        missing.append("evolving_geodesic.json")
    if not c.get("emit_sweep"):
        missing.append("emit_sweep (>=2 launches)")
    if not c.get("psi4_l2m0"):
        missing.append("psi4_mode_l2m0.dat with data rows")
    if not c.get("psi4_l2_all"):
        missing.append("psi4_mode_l2_all.dat")
    if not c.get("score"):
        missing.append("score.json")
    if summary.get("f_geo") is None:
        missing.append("f_geo")
    if summary.get("h_quality_ok") is not True:
        missing.append("h_quality_ok=true")
    return missing


def analyze_resolution_ladder(
    runs: dict[str, Path],
    spacings: dict[str, float],
    *,
    rel_tol: float = 0.10,
) -> dict[str, Any]:
    """Compare coarse/medium/fine summaries and Richardson-extrapolate f_geo."""
    summaries = {k: extract_run_summary(p) for k, p in runs.items()}
    f_geo = {k: float(summaries[k]["f_geo"]) for k in runs if summaries[k].get("f_geo") is not None}
    report: dict[str, Any] = {
        "summaries": summaries,
        "spacings": spacings,
        "f_geo": f_geo,
        "richardson_f_geo": None,
        "fine_vs_medium_rel": None,
        "pass_rel_tol": None,
        "missing": {k: assess_completeness_for_validation(summaries[k]) for k in runs},
    }
    if set(f_geo) >= {"coarse", "medium", "fine"}:
        report["richardson_f_geo"] = unequal_richardson(f_geo, spacings)
        rel = relative_diff(f_geo["fine"], f_geo["medium"])
        report["fine_vs_medium_rel"] = rel
        report["pass_rel_tol"] = rel is not None and rel <= rel_tol
    return report


def analyze_domain_ladder(
    runs: dict[str, Path],
    *,
    reference: str = "medium",
    large: str = "large",
    rel_tol: float = 0.10,
) -> dict[str, Any]:
    summaries = {k: extract_run_summary(p) for k, p in runs.items()}
    report: dict[str, Any] = {
        "summaries": summaries,
        "missing": {k: assess_completeness_for_validation(summaries[k]) for k in runs},
        "ref_vs_large_rel": None,
        "pass_rel_tol": None,
    }
    fa = summaries.get(reference, {}).get("f_geo")
    fb = summaries.get(large, {}).get("f_geo")
    if fa is not None and fb is not None:
        rel = relative_diff(float(fa), float(fb))
        report["ref_vs_large_rel"] = rel
        report["pass_rel_tol"] = rel is not None and rel <= rel_tol
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_one = sub.add_parser("summarize", help="Summarize one HQ run")
    p_one.add_argument("run_dir", type=Path)
    p_one.add_argument("--out", type=Path, default=None)

    p_check = sub.add_parser(
        "check-complete",
        help="Exit 1 if run is incomplete for emission-sweep/GW validation",
    )
    p_check.add_argument("run_dir", type=Path)

    p_res = sub.add_parser("resolution", help="Richardson ladder for 3 runs")
    p_res.add_argument("--coarse", type=Path, required=True)
    p_res.add_argument("--medium", type=Path, required=True)
    p_res.add_argument("--fine", type=Path, required=True)
    p_res.add_argument("--h-coarse", type=float, required=True)
    p_res.add_argument("--h-medium", type=float, required=True)
    p_res.add_argument("--h-fine", type=float, required=True)
    p_res.add_argument("--rel-tol", type=float, default=0.10)
    p_res.add_argument("--out", type=Path, default=None)

    p_dom = sub.add_parser("domain", help="Domain-size comparison")
    p_dom.add_argument("--small", type=Path, default=None)
    p_dom.add_argument("--medium", type=Path, required=True)
    p_dom.add_argument("--large", type=Path, required=True)
    p_dom.add_argument("--rel-tol", type=float, default=0.10)
    p_dom.add_argument("--out", type=Path, default=None)

    args = parser.parse_args()

    if args.cmd == "summarize":
        summary = extract_run_summary(args.run_dir)
        summary["validation_missing"] = assess_completeness_for_validation(summary)
        text = json.dumps(summary, indent=2)
        if args.out:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(text + "\n", encoding="utf-8")
        print(text)
        return 0

    if args.cmd == "check-complete":
        summary = extract_run_summary(args.run_dir)
        missing = assess_completeness_for_validation(summary)
        payload = {"path": str(args.run_dir), "missing": missing, "ok": not missing}
        print(json.dumps(payload, indent=2))
        return 0 if not missing else 1

    if args.cmd == "resolution":
        report = analyze_resolution_ladder(
            {
                "coarse": args.coarse,
                "medium": args.medium,
                "fine": args.fine,
            },
            {
                "coarse": args.h_coarse,
                "medium": args.h_medium,
                "fine": args.h_fine,
            },
            rel_tol=args.rel_tol,
        )
        text = json.dumps(report, indent=2)
        if args.out:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(text + "\n", encoding="utf-8")
        print(text)
        return 0

    if args.cmd == "domain":
        runs = {"medium": args.medium, "large": args.large}
        if args.small is not None:
            runs["small"] = args.small
        report = analyze_domain_ladder(runs, rel_tol=args.rel_tol)
        text = json.dumps(report, indent=2)
        if args.out:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(text + "\n", encoding="utf-8")
        print(text)
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
