#!/usr/bin/env python3
"""Derive the article tables from the packed per-cell time series.

Reads ``data/<cell>/{sector_barycenters,confinement}.dat`` + ``metadata.json``
and writes:

  analysis/trajectories.csv  -- drift/separation/core series, every cell
  analysis/summary.csv       -- one row per cell (birth checks -> final state)
  analysis/summary.md        -- the same table, ready to paste into the article

Pure stdlib: the packed streams are a few hundred rows each.

Usage: python3 analysis/make_tables.py [PACK_DIR]   (default: parent of this file)
"""

from __future__ import annotations

import csv
import json
import math
import pathlib
import sys

# Column indices (0-based) in the GRTeclyn small_data streams.
SB = {"t": 0, "tot_c": 1, "x_c": 2, "rms_c": 5, "tot_p": 6, "x_p": 7, "rms_p": 10}
CF = {"t": 0, "peak": 2, "conf": 4, "min_chi": 17}

# cell -> (sectors present, geometry description)
# The dx = 0.5 production runs that the convergence campaign never twinned:
# the single stars and the mirrored pair, plus the unmirrored pair they are
# compared against.  Everything else that used to live here was superseded by
# the convA_* ladder and is covered by convergence_check.py instead.
CELLS = [
    ("single_p", "canonical", "1 star at grid centre"),
    ("single_m", "phantom", "1 star at grid centre"),
    ("pair_pm", "canonical + phantom", "x = +4 / -4, sep 8"),
    ("pair_mp_mirror", "phantom + canonical (mirror control)", "x = +4 / -4, sep 8, sectors swapped"),
]
SAMPLE_DT = 4.0  # trajectories.csv row spacing in code time


def read_stream(path: pathlib.Path) -> list[list[float]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line[0].isalpha():
            continue
        rows.append([float(v) for v in line.split()])
    return rows


def at_time(rows: list[list[float]], t: float) -> list[float] | None:
    """Nearest sample to time t (streams are dense and uniformly spaced)."""
    if not rows:
        return None
    return min(rows, key=lambda r: abs(r[0] - t))


def finite(value: float) -> float | str:
    return "" if value is None or math.isnan(value) else value


def main(argv: list[str]) -> int:
    pack = pathlib.Path(argv[0]) if argv else pathlib.Path(__file__).resolve().parent.parent
    data = pack / "campaign"
    out = pack / "analysis"
    out.mkdir(parents=True, exist_ok=True)

    traj_rows: list[dict] = []
    summary_rows: list[dict] = []

    for cell, sectors, geometry in CELLS:
        cdir = data / cell
        if not cdir.is_dir():
            continue
        sb = read_stream(cdir / "sector_barycenters.dat")
        cf = read_stream(cdir / "confinement.dat")
        meta = json.loads((cdir / "metadata.json").read_text(encoding="utf-8"))
        conv = meta.get("grtresna_convergence", {})
        if not sb:
            continue

        first, last = sb[0], sb[-1]
        t_end = last[SB["t"]]

        # --- trajectory series, sampled every SAMPLE_DT ---
        t = 0.0
        while t <= t_end + 1e-9:
            row = at_time(sb, t)
            core = at_time(cf, t)
            x_c, x_p = row[SB["x_c"]], row[SB["x_p"]]
            sep = x_c - x_p if not (math.isnan(x_c) or math.isnan(x_p)) else float("nan")
            traj_rows.append({
                "cell": cell,
                "t": round(row[SB["t"]], 3),
                "x_canon": finite(x_c),
                "x_phantom": finite(x_p),
                "separation": finite(sep),
                "drift_canon": finite(x_c - first[SB["x_c"]]),
                "drift_phantom": finite(x_p - first[SB["x_p"]]),
                "rms_canon": finite(row[SB["rms_c"]]),
                "rms_phantom": finite(row[SB["rms_p"]]),
                "peak_amplitude": finite(core[CF["peak"]]) if core else "",
                "confined_frac": finite(core[CF["conf"]]) if core else "",
                "min_chi": finite(core[CF["min_chi"]]) if core else "",
            })
            t += SAMPLE_DT

        # --- per-cell summary ---
        core0, core1 = at_time(cf, 0.0), (cf[-1] if cf else None)
        drift_c = last[SB["x_c"]] - first[SB["x_c"]]
        drift_p = last[SB["x_p"]] - first[SB["x_p"]]
        sep0 = first[SB["x_c"]] - first[SB["x_p"]]
        sep1 = last[SB["x_c"]] - last[SB["x_p"]]
        summary_rows.append({
            "cell": cell,
            "sectors": sectors,
            "geometry": geometry,
            "t_end": round(t_end, 2),
            "ham_pct": conv.get("ham_pct", ""),
            "mom_pct": conv.get("mom_pct", ""),
            "t0_total_canon": finite(first[SB["tot_c"]]),
            "t0_rms_canon": finite(first[SB["rms_c"]]),
            "t0_total_phantom": finite(first[SB["tot_p"]]),
            "t0_rms_phantom": finite(first[SB["rms_p"]]),
            "drift_canon": finite(drift_c),
            "drift_phantom": finite(drift_p),
            "sep_initial": finite(sep0),
            "sep_final": finite(sep1),
            "peak_t0": finite(core0[CF["peak"]]) if core0 else "",
            "peak_end": finite(core1[CF["peak"]]) if core1 else "",
            "min_chi_t0": finite(core0[CF["min_chi"]]) if core0 else "",
            "min_chi_end": finite(core1[CF["min_chi"]]) if core1 else "",
        })

    with (out / "trajectories.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(traj_rows[0].keys()))
        writer.writeheader()
        writer.writerows(traj_rows)

    with (out / "summary.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    # Markdown mirror of summary.csv for direct paste into the article.
    def fmt(value, digits=3):
        return f"{value:.{digits}f}" if isinstance(value, float) else (value or "—")

    lines = [
        "# Per-cell summary (generated by analysis/make_tables.py)",
        "",
        "| cell | sectors | t_end | Ham % | Mom % | drift canon | drift phantom | sep 0 -> end | min chi 0 -> end |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in summary_rows:
        sep = (
            f"{fmt(row['sep_initial'], 2)} -> {fmt(row['sep_final'], 2)}"
            if row["sep_initial"] != "" else "—"
        )
        lines.append(
            f"| `{row['cell']}` | {row['sectors']} | {row['t_end']:.0f} | "
            f"{fmt(row['ham_pct'], 4)} | {fmt(row['mom_pct'], 4)} | "
            f"{fmt(row['drift_canon'])} | {fmt(row['drift_phantom'])} | {sep} | "
            f"{fmt(row['min_chi_t0'], 4)} -> {fmt(row['min_chi_end'], 4)} |"
        )
    (out / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[tables] trajectories.csv: {len(traj_rows)} rows")
    print(f"[tables] summary.csv: {len(summary_rows)} cells")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
