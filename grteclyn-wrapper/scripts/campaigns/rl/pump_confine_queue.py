#!/usr/bin/env python3
"""Confinement campaign A: can the pump HOLD matter instead of merely slowing it?

Six arms, all cloned from the same RM baseline as runs/pump_ladder_m0 and all
with the pump always on (rl_pump_stop_time = -1), so every arm is directly
comparable to lad_m0_tp30. Only the knob named in each arm differs.

WHAT THE LADDER MEASURED (runs/pump_ladder_m0, confined_frac by sector):

  * canonical -- the pump does not leak it, the pump EJECTS it at startup.
    tp30 loses 0.143 (21%) of confined_frac_canon in t < 5.8 while the
    pump-free rung loses 0.029; after t ~ 13 the pumped run is flat to three
    decimals (lambda ~ 0.000) and even recovers to 0.504 by t = 30. The whole
    canonical deficit is a startup transient, not a slow leak.
  * phantom -- no startup transient at all (0.7841 pumped vs 0.7836 free at
    t = 5.76), but a decay that is only slowed: lambda 0.217 free -> 0.076
    pumped over [13, 30]. Never arrested.

Two different failures, so two different arms. Gains and grip range attack the
phantom decay; the canonical well_depth bracket attacks the startup ejection
(the target amplitude is a flat 0.15 on all five lumps, canonical and phantom
alike, which is the prime suspect for the sector asymmetry).

The clean measurement window is [13, 27]: tp30's safety governor holds at
1.0000 until t ~ 28.9, then a late L2_Ham runaway (6.0e-3 at t=27.5 -> 4.8e-2
at t=30) chokes it to 2.2e-4. Any arm whose governor crosses 0.5 before t = 27
has bought its confinement with constraint violation and is disqualified.

Requires a binary with the per-sector gains (commit 9e81d480 or later); check
with `strings main3d.gnu.CUDA.ex | grep rl_pump_kp_phantom`.

Inherits the mandatory node-local plotfile scratch pattern from
pump_ladder_queue.py -- see scripts/campaigns/README.md and the wrapper README
section "Plotfile scratch MUST be node-local". Deletion is ledger-gated, so
every extraction this campaign will ever need must be enabled at launch:
confinement (the pass/fail metric) and the metric_stack cache for the post-run
4D geodesic score. Nothing can be recovered afterwards.
"""

from __future__ import annotations

import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/home/jovyan/nachevsky/test/simulation/GRTeclyn")
BASE = ROOT / "runs/grtresna_promote/bcma_rm_L128_N256_t30_hq_eval000146/params.txt"
# Overridable with --out; campaign B (corrected spotlight routing) writes to
# runs/pump_confine_b so the broken-routing A record stays separate.
OUT = ROOT / "runs/pump_confine_a"
EXE = ROOT / "Examples/RadialRecipe/main3d.gnu.CUDA.ex"
WRAPPER = ROOT / "grteclyn-wrapper"
# Call the project venv python DIRECTLY rather than via `uv run`, which would
# write to ~/.cache/uv.
VENV_PY = WRAPPER / ".venv/bin/python"
SCORE_EVOLVING = WRAPPER / "scripts/campaigns/rl/score_evolving_geodesic.py"

# Heavy transients stay on node-local NVMe; only .dat + small_data reach NFS.
SCRATCH = Path("/tmp/grteclyn_scratch")

# The wrapper README documents 6 concurrent HQ runs, which was the GPU count on
# the ladder node, not an I/O ceiling -- with node-local scratch the limit is
# GPUs. This node has 8 (plus 192 cores / 1.4 TB RAM for the 8 sidecar
# consumers), so all 8 arms run in a single batch.
MAX_CONCURRENT = 8
GPUS = [0, 1, 2, 3, 4, 5, 6, 7]
PLOT_INTERVAL = 144
DRAIN_SECONDS = 600
KEEP_LAST = 3

# Steady state is n_runs * (keep_last + jobs) * plotfile_size
# = 8 * 4 * ~2.9 GB ~= 93 GB. Refuse to start without headroom over that.
MIN_FREE_GB = 140

# Canonical spotlights are lumps 0 and 3; phantom are 1, 2, 4
# (recipe_scalar_field_signs = 1 -1 -1 1 -1).
#
# (name, [param lines])  -- every key listed here is stripped from the baseline
# before being re-appended, so each appears exactly once in the final file.
ARMS = [
    # -- phantom authority ladder. Prediction: lambda_phan over [13,27] drops
    #    below the 0.075 baseline, WITHOUT moving confined_frac_canon.
    #    Stability: omega_eff^2 ~ k_p*w, so k_p=12 gives omega*dt ~ 0.055 at the
    #    finest dt ~ 0.016; x4 gives ~0.11. Do not exceed x4 without lowering
    #    dt_multiplier.
    ("pca_pg2", [
        "rl_pump_kp_phantom = 24.0",
        "rl_pump_kd_phantom = 14.0",
    ]),
    ("pca_pg4", [
        "rl_pump_kp_phantom = 48.0",
        "rl_pump_kd_phantom = 28.0",
    ]),
    # -- grip range. The PD envelope is sech(r/1.667), which is ~4e-5 at the
    #    escaped-matter radius r ~ 18: gain alone cannot reach matter that has
    #    already left the well. Watch rms_radius_* -- the envelope also scales
    #    the TARGET amplitude, so a wider held blob is not a win.
    ("pca_tw2", [
        "rl_pump_target_width = 3.333333333333334",
    ]),
    # -- do authority and grip compose, or is one of them the whole story?
    ("pca_pg4tw2", [
        "rl_pump_kp_phantom = 48.0",
        "rl_pump_kd_phantom = 28.0",
        "rl_pump_target_width = 3.333333333333334",
    ]),
    # -- canonical startup transient, bracketed. The target amplitude is a flat
    #    well_depth = 0.15 everywhere; if the canonical initial state sits above
    #    it the PD drags the field down and ejects matter, if below it kicks.
    #    Whichever direction shrinks the t<5.8 loss identifies the sign.
    ("pca_cwd_up", [
        "trajectory_lump0_well_depth = 0.30",
        "trajectory_lump3_well_depth = 0.30",
    ]),
    ("pca_cwd_dn", [
        "trajectory_lump0_well_depth = 0.075",
        "trajectory_lump3_well_depth = 0.075",
    ]),
    # -- the same bracket on the phantom side. The ladder shows phantom has no
    #    startup transient, which is the evidence that 0.15 is roughly right
    #    for that sector -- but that is an inference, not a measurement, and
    #    target mismatch is a different lever from gain (arms 1-2) and grip
    #    (arm 3). Symmetric brackets make the sector comparison clean.
    ("pca_pwd_up", [
        "trajectory_lump1_well_depth = 0.30",
        "trajectory_lump2_well_depth = 0.30",
        "trajectory_lump4_well_depth = 0.30",
    ]),
    ("pca_pwd_dn", [
        "trajectory_lump1_well_depth = 0.075",
        "trajectory_lump2_well_depth = 0.075",
        "trajectory_lump4_well_depth = 0.075",
    ]),
    # -- SECOND WAVE (launched after the two _up arms were retired; see below).
    #    The seed lumps have central amplitude phi_c = 0.08
    #    (initial_data.matter.json), while the pump commands well_depth = 0.15
    #    -- it drives every lump to ~1.9x its true height. These two arms set
    #    the target to the matter that actually exists, on BOTH sectors at once
    #    (the _dn arms each fix one sector and leave the other over-driven).
    ("pca_match", [
        "trajectory_lump0_well_depth = 0.08",
        "trajectory_lump1_well_depth = 0.08",
        "trajectory_lump2_well_depth = 0.08",
        "trajectory_lump3_well_depth = 0.08",
        "trajectory_lump4_well_depth = 0.08",
    ]),
    # Matched target PLUS maximum phantom authority. The open question the
    # ladder leaves is why the pump injects 2.7x into the canonical sector and
    # exactly 0x into the phantom one despite identical target, gains and seed
    # amplitude. If the phantom sector is authority-limited this arm moves it;
    # if it is limited by its own self-repulsion, nothing here will.
    ("pca_match_pg4", [
        "trajectory_lump0_well_depth = 0.08",
        "trajectory_lump1_well_depth = 0.08",
        "trajectory_lump2_well_depth = 0.08",
        "trajectory_lump3_well_depth = 0.08",
        "trajectory_lump4_well_depth = 0.08",
        "rl_pump_kp_phantom = 48.0",
        "rl_pump_kd_phantom = 28.0",
    ]),
]

# RETIRED 2026-07-28, kept here as the record of a falsified direction:
# pca_cwd_up / pca_pwd_up raised well_depth 0.15 -> 0.30. Both refined ~111x
# more level-3 cells than their _dn counterparts (12.3M vs 110k), ran ~3x
# slower, and pca_pwd_up aborted with "NaN in K" at t=4.31. Over-driving the
# target does not merely waste grid, it destroys the evolution.
#
# CAMPAIGN A POSTMORTEM (2026-07-28, Debug.md 19.8): every pca_* arm above ran
# on a binary whose recipe_scalar_field_signs parser kept only the FIRST value
# of "1 -1 -1 1 -1" -- all five spotlights drove the canonical field, none
# drove Phi-. pca_pg2/pg4 came out bit-identical to lad_m0_tp30 (the phantom
# gains multiplied a force that never existed), which is how the bug was
# caught. Campaign B below reruns the decisive arms on the fixed binary; its
# first-slice check is the "recipe_scalar_field_signs parsed: 1 -1 -1 1 -1"
# echo in run.log.
ARMS_B = [
    # Corrected baseline: stock tp30 knobs, routing fixed -- the first run in
    # which the pump ACTUALLY drives the phantom sector (2 canonical sites +
    # 3 phantom sites, as configured all along).
    ("pcb_base", []),
    # Phantom authority x4 on top of correct routing -- the original 18.6
    # question, now actually testable.
    ("pcb_pg4", [
        "rl_pump_kp_phantom = 48.0",
        "rl_pump_kd_phantom = 28.0",
    ]),
    # Target matched to the seed lumps' true amplitude (phi_c = 0.08).
    ("pcb_match", [
        "trajectory_lump0_well_depth = 0.08",
        "trajectory_lump1_well_depth = 0.08",
        "trajectory_lump2_well_depth = 0.08",
        "trajectory_lump3_well_depth = 0.08",
        "trajectory_lump4_well_depth = 0.08",
    ]),
    ("pcb_match_pg4", [
        "trajectory_lump0_well_depth = 0.08",
        "trajectory_lump1_well_depth = 0.08",
        "trajectory_lump2_well_depth = 0.08",
        "trajectory_lump3_well_depth = 0.08",
        "trajectory_lump4_well_depth = 0.08",
        "rl_pump_kp_phantom = 48.0",
        "rl_pump_kd_phantom = 28.0",
    ]),
]
ARMS += ARMS_B

# --- campaign C -------------------------------------------------------------
# Campaign B's verdict, read off the first ~19 time units: the TARGET AMPLITUDE
# is the whole story and phantom GAIN is not a lever. pcb_base and pcb_pg4 (x4
# phantom authority) trace the SAME collapse curve to three digits -- min_chi
# 0.133 vs 0.139 at t=8.64 -- and pcb_base died at t=12.27 with "NaN in K" after
# min_chi fell to 1.6e-3, the same signature as campaign A's pwd_up. Both
# over-target arms inject 48 -> ~88 absolute in the first 3 units and then
# collapse under the mass they added. The matched arms (target 0.08 = the seed
# lumps' true phi_c) sit flat at ~44 with min_chi 0.64 at t=18.7.
#
# So campaign C brackets the surviving question: how much injection is safe, and
# does grip range recover the ~19% phantom shed in the first 1.4 units?
ARMS_C = [
    # Injection ceiling. 0.08 is stable, 0.15 is lethal -- bisect once. If this
    # survives to t=30 the safe target is >= 0.10 and there is headroom to gain
    # matter; if it collapses, matched-only is the operating point.
    ("pcc_t010", [
        "trajectory_lump0_well_depth = 0.10",
        "trajectory_lump1_well_depth = 0.10",
        "trajectory_lump2_well_depth = 0.10",
        "trajectory_lump3_well_depth = 0.10",
        "trajectory_lump4_well_depth = 0.10",
    ]),
    # Grip range x2 at the safe target. The sech(r/w) envelope is ~4e-5 at the
    # escaped-matter radius, so widening is the only lever that can reach matter
    # already outside the well. Watch rms_radius_* (cols 13/15): a wider held
    # blob is not a win.
    ("pcc_match_tw2", [
        "trajectory_lump0_well_depth = 0.08",
        "trajectory_lump1_well_depth = 0.08",
        "trajectory_lump2_well_depth = 0.08",
        "trajectory_lump3_well_depth = 0.08",
        "trajectory_lump4_well_depth = 0.08",
        "rl_pump_target_width = 3.333333333333334",
    ]),
]
ARMS += ARMS_C

# --- campaign D -------------------------------------------------------------
# Two arms launched onto the GPUs freed by the over-target deaths (pcb_base NaN
# at t=12.27; pcb_pg4 killed at t=10.08 on the min_chi < 0.05 criterion, having
# tracked pcb_base to three digits the whole way down).
#
# The gains here are honoured as LITERAL ZERO, not "inherit": RLPumpForce.hpp:62
# selects the phantom gain on `k_p_phantom >= 0.0`, and only a NEGATIVE value
# falls back to the canonical k_p. So 0.0 really does mean "no PD force on the
# phantom sites" -- which is the control this project has never had.
ARMS_D = [
    # THE DECISIVE CONTROL. Matched target, canonical sites pumped normally,
    # phantom sites pumped with zero gain. Sec. 19.8 concluded that every
    # "phantom benefits from the pump" result to date was INDIRECT (geometric
    # shielding by over-pumped canonical matter), because the phantom sites were
    # mis-routed and received no force at all. This arm reproduces that
    # condition deliberately at a survivable target: if the phantom sector still
    # holds flat here, the direct forcing in pcb_match is not what confines it,
    # and the confinement claim belongs to geometry, not the controller.
    ("pcd_match_p0", [
        "trajectory_lump0_well_depth = 0.08",
        "trajectory_lump1_well_depth = 0.08",
        "trajectory_lump2_well_depth = 0.08",
        "trajectory_lump3_well_depth = 0.08",
        "trajectory_lump4_well_depth = 0.08",
        "rl_pump_kp_phantom = 0.0",
        "rl_pump_kd_phantom = 0.0",
    ]),
    # Persistence. "Confined" has only ever been measured to t=30, and the tp30
    # ladder rung choked on a constraint runaway at t~29 -- so t=30 cannot
    # distinguish "held" from "not yet dispersed". Double the horizon at the
    # matched target and find out which.
    ("pcd_match_t60", [
        "trajectory_lump0_well_depth = 0.08",
        "trajectory_lump1_well_depth = 0.08",
        "trajectory_lump2_well_depth = 0.08",
        "trajectory_lump3_well_depth = 0.08",
        "trajectory_lump4_well_depth = 0.08",
        "stop_time = 60.0",
    ]),
]
ARMS += ARMS_D

# --- campaign E -------------------------------------------------------------
# The superposed-target law (rl_pump_superpose_targets, RLPumpForce.hpp).
# Diagnosis from campaign B geometry: the legacy per-site PD errors STRIP the
# overlap between same-sector sites -- each site treats its neighbour's lump as
# excess and their down-drives sum.  Canonical sites sit 8.9 apart (overlap
# 9e-3, sector +11% by t=1.44); the three phantom sites sit 4.6-5.0 apart
# (overlap ~0.1, sector -19%).  With the law ON, a sector's sites share one
# summed target and a capped weight; identical to legacy for isolated sites.
# A/B twins: pce_sup vs pcb_match (t=30), pce_sup_t60 vs pcd_match_t60 (t=60).
ARMS_E = [
    ("pce_sup", [
        "trajectory_lump0_well_depth = 0.08",
        "trajectory_lump1_well_depth = 0.08",
        "trajectory_lump2_well_depth = 0.08",
        "trajectory_lump3_well_depth = 0.08",
        "trajectory_lump4_well_depth = 0.08",
        "rl_pump_superpose_targets = 1",
    ]),
    ("pce_sup_t60", [
        "trajectory_lump0_well_depth = 0.08",
        "trajectory_lump1_well_depth = 0.08",
        "trajectory_lump2_well_depth = 0.08",
        "trajectory_lump3_well_depth = 0.08",
        "trajectory_lump4_well_depth = 0.08",
        "rl_pump_superpose_targets = 1",
        "stop_time = 60.0",
    ]),
]
ARMS += ARMS_E

# --- campaign F -------------------------------------------------------------
# Campaign C/D/E verdict (Debug.md 19.13): pcc_t010 (aim 0.10) is the best
# configuration found -- healthiest geometry by 6x in min_chi, exotic sector
# essentially conserved (-5%), and no early strip at all, because aiming
# slightly ABOVE the matter's own amplitude stops the legacy per-site law
# reading a neighbour's lump as excess. But pcd_match_t60 then showed that
# pcb_match's "clean finish" at t=30 was collapse already in progress: extended
# past its twin's stop time it died at t~32. So t=30 proves nothing and the
# standing evidence rule is now t >= 40.
#
# Every arm here therefore runs to t=60, and the watchdog (see below) aborts a
# collapsing arm early so the GPU is not spent on a dead run. All four clone
# pcc_t010 unless noted. recipe_scalar_field_signs = 1 -1 -1 1 -1, so lumps
# 0 and 3 are canonical and 1, 2, 4 are exotic -- and well_depth is already per
# lump, so per-sector aiming needs no code change.
ARMS_F = [
    # PRIMARY. Does the best-known config survive past t~32, where aim 0.08
    # died? This is also the field validation of the theta_plus refinement-edge
    # fix: pcc_t010 on the old binary reported theta<0 from t=25.64 with the
    # minimum pinned at a fixed r=10.40 while its interior stayed healthy
    # (lapse 0.601, chi 0.499 at t=30).
    ("pcf_t010_t60", [
        "trajectory_lump0_well_depth = 0.10",
        "trajectory_lump1_well_depth = 0.10",
        "trajectory_lump2_well_depth = 0.10",
        "trajectory_lump3_well_depth = 0.10",
        "trajectory_lump4_well_depth = 0.10",
        "stop_time = 60.0",
    ]),
    # Is the stability optimum above 0.10, or is the 0.15 cliff already near?
    # 0.15 injects monotonically and NaNs at t~12; 0.10 conserves. Bisect up.
    ("pcf_t012_t60", [
        "trajectory_lump0_well_depth = 0.12",
        "trajectory_lump1_well_depth = 0.12",
        "trajectory_lump2_well_depth = 0.12",
        "trajectory_lump3_well_depth = 0.12",
        "trajectory_lump4_well_depth = 0.12",
        "stop_time = 60.0",
    ]),
    # Per-sector aiming. Exotic matter is what resists collapse (pcd_match_p0
    # lost its NEC violation entirely and then collapsed hardest of any arm),
    # and aiming high is what conserves a sector -- but aiming high on the
    # canonical sector is what over-fills it. So aim high only where it pays.
    ("pcf_split_t60", [
        "trajectory_lump0_well_depth = 0.08",
        "trajectory_lump3_well_depth = 0.08",
        "trajectory_lump1_well_depth = 0.12",
        "trajectory_lump2_well_depth = 0.12",
        "trajectory_lump4_well_depth = 0.12",
        "stop_time = 60.0",
    ]),
    # The superposed law (19.11) eliminated the overlap strip and then killed
    # its run by over-feeding the bridges between the exotic chain (67 units of
    # mass instead of 44, collapse at t~14-19). Reduced aim is the test of
    # whether that was the law or the dose. Pass = total ~48, not ~67, at t=13,
    # with the exotic strip still absent at t=1.44.
    ("pcf_sup_a06_t60", [
        "trajectory_lump0_well_depth = 0.06",
        "trajectory_lump1_well_depth = 0.06",
        "trajectory_lump2_well_depth = 0.06",
        "trajectory_lump3_well_depth = 0.06",
        "trajectory_lump4_well_depth = 0.06",
        "rl_pump_superpose_targets = 1",
        "stop_time = 60.0",
    ]),
]
ARMS += ARMS_F

# Stripped from the baseline for every arm. amr.plot_file / amr.check_file are
# absolute paths INTO THE BASELINE RUN DIR: if they survive the clone, this
# campaign writes plotfiles into the baseline and its --delete consumer prunes
# the baseline's own data.
STRIP_ALWAYS = {
    "output_path", "plot_interval", "rl_pump_stop_time",
    "controller_reservoir_mode", "amr.plot_file", "amr.check_file",
}


# --- collapse watchdog ------------------------------------------------------
# Kill criteria (Debug.md 19.6, extended in 19.14). min_lapse is the new
# PRIMARY trigger: in both confirmed collapses (pcd_match_t60, pce_sup) it led
# min_chi by ~1 time unit, so the old min_chi<0.05 line is a LATE indicator that
# spends a GPU on a run already dead. Constraint criteria only apply before
# t=27, where the tp30 baseline's own governor is still 1.0000; after that a
# choke is expected and is not evidence about the arm.
WATCH_MIN_LAPSE = 0.15
WATCH_MIN_CHI = 0.05
WATCH_L2_HAM = 0.035
WATCH_GOVERNOR = 0.5
WATCH_CONSTRAINT_UNTIL_T = 27.0

# Column layout is positional -- these .dat files carry no header row. Indices
# are 0-based and must track the writers in Examples/RadialRecipe/
# RadialRecipeLevel.cpp (write_header_line calls). N_COLS is checked on every
# read so a column-layout change fails the watchdog OPEN (no kill) instead of
# silently reading the wrong quantity and aborting healthy runs.
COLLAPSE_N_COLS = 15   # time min_lapse min_chi max_abs_K ... pump_work
COLLAPSE_I_TIME, COLLAPSE_I_LAPSE, COLLAPSE_I_CHI = 0, 1, 2
NORMS_N_COLS = 11      # time L2_Ham L2_Mom ... pump_force_L2 governor pump_fi_L2
NORMS_I_TIME, NORMS_I_HAM, NORMS_I_GOV = 0, 1, 9


def last_row(path: Path, n_cols: int) -> list[float] | None:
    """Last complete row of a positional .dat file, or None if unreadable.

    Returns None on a short/absent/partially-written file: the sim appends
    while we poll, so a torn final line is normal and must not be a kill.
    """
    try:
        lines = path.read_text().splitlines()
    except OSError:
        return None
    for ln in reversed(lines):
        fields = ln.split()
        if len(fields) != n_cols:
            continue
        try:
            return [float(f) for f in fields]
        except ValueError:
            continue
    return None


def collapse_verdict(d: Path) -> str | None:
    """Reason to abort this run now, or None to let it continue."""
    row = last_row(d / "data/collapse_diagnostics.dat", COLLAPSE_N_COLS)
    if row is not None:
        t = row[COLLAPSE_I_TIME]
        lapse, chi = row[COLLAPSE_I_LAPSE], row[COLLAPSE_I_CHI]
        if any(v != v or v in (float("inf"), float("-inf")) for v in row):
            return f"non-finite value in collapse_diagnostics at t={t:.2f}"
        if lapse < WATCH_MIN_LAPSE:
            return (f"min_lapse={lapse:.4g} < {WATCH_MIN_LAPSE} at t={t:.2f} "
                    f"(min_chi={chi:.4g})")
        if chi < WATCH_MIN_CHI:
            return f"min_chi={chi:.4g} < {WATCH_MIN_CHI} at t={t:.2f}"

    row = last_row(d / "data/constraint_norms.dat", NORMS_N_COLS)
    if row is not None:
        t = row[NORMS_I_TIME]
        ham, gov = row[NORMS_I_HAM], row[NORMS_I_GOV]
        if ham != ham:
            return f"L2_Ham is NaN at t={t:.2f}"
        if t < WATCH_CONSTRAINT_UNTIL_T:
            if ham > WATCH_L2_HAM:
                return f"L2_Ham={ham:.4g} > {WATCH_L2_HAM} at t={t:.2f}"
            if gov < WATCH_GOVERNOR:
                return f"governor={gov:.4g} < {WATCH_GOVERNOR} at t={t:.2f}"
    return None


def param_key(line: str) -> str | None:
    """Key of a params.txt assignment, or None for blanks/comments."""
    s = line.strip()
    if not s or s.startswith("#") or "=" not in s:
        return None
    return s.split("=", 1)[0].strip()


def log(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "queue.log").open("a") as fh:
        fh.write(line + "\n")


def scratch_dir(name: str) -> Path:
    return SCRATCH / name


def prepare(name: str, extra: list[str]) -> Path:
    d = OUT / name
    (d / "data").mkdir(parents=True, exist_ok=True)
    scratch = scratch_dir(name)
    scratch.mkdir(parents=True, exist_ok=True)

    strip = set(STRIP_ALWAYS)
    for line in extra:
        key = param_key(line)
        if key:
            strip.add(key)

    keep = [line for line in BASE.read_text().splitlines()
            if param_key(line) not in strip]
    keep += [
        "",
        f"# --- confinement campaign A: {name} ---",
        f'output_path = "{d}"',
        f'amr.check_file = "{scratch}/RadialRecipeChk"',
        f'amr.plot_file  = "{scratch}/RadialRecipePlt"',
        f"plot_interval       = {PLOT_INTERVAL}",
        "rl_pump_stop_time = -1.0",
        "controller_reservoir_mode = 0",
        *extra,
    ]
    (d / "params.txt").write_text("\n".join(keep) + "\n")
    return d


def start_consumer(d: Path) -> subprocess.Popen:
    # Deletion is gated on the extraction ledger, so anything not enabled here
    # is unrecoverable once scratch is purged. --confinement-timeseries is the
    # campaign's pass/fail metric; --evolving-geodesic builds the metric_stack
    # cache the post-run scorer needs.
    cmd = [
        str(VENV_PY), "-m",
        "grteclyn_wrapper.visualisation.process_wave.consume_plotfiles",
        "--data", str(scratch_dir(d.name)), "--out", str(d / "small_data"),
        "--center", "64", "64", "64",
        "--radii", "8", "12", "24", "--n-points", "32",
        "--no-psi4", "--ftl-timeseries", "--ftl-l", "8",
        "--evolving-geodesic",
        "--confinement-timeseries", "--confinement-well-width", "1.667",
        "--watch", "--delete", "--keep-last", str(KEEP_LAST),
        "--stable-seconds", "30", "--poll-seconds", "20", "-j", "1", "--verbose",
    ]
    fh = (d / "consumer.log").open("w")
    return subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT,
                            env=consumer_env(), start_new_session=True)


def consumer_env() -> dict:
    """Environment shared by the consumer AND the post-run scoring pass.

    They must agree: the scorer reads the cache the consumer wrote, so a
    mismatch in GRTECLYN_METRIC_STACK_N_SPACE would validate against the wrong
    resolution. 257 = 2*half_width/finest_dx + 1, matching the run's finest AMR
    level; at the default 33 the cache resamples at the COARSEST dx and erases
    everything on the refined levels.
    """
    env = dict(os.environ)
    env["GEODESIC_EMIT_MIN_TIME"] = "0"
    env.pop("RL_PUMP_STOP_TIME", None)
    env.pop("GRTECLYN_KEEP_PLOTFILES", None)  # must not disable deletion
    env["GRTECLYN_EVOLVING_GEODESIC_MODE"] = "hq"
    env["GRTECLYN_METRIC_STACK_N_SPACE"] = "257"
    # Keep every cache/temp write inside our own scratch, never ~/.local or
    # ~/.cache (~/.local is admin-locked on this node).
    cache = SCRATCH / "_cache"
    cache.mkdir(parents=True, exist_ok=True)
    env["UV_CACHE_DIR"] = str(cache / "uv")
    env["XDG_CACHE_HOME"] = str(cache)
    env["MPLCONFIGDIR"] = str(cache / "mpl")
    env["TMPDIR"] = str(cache / "tmp")
    (cache / "tmp").mkdir(parents=True, exist_ok=True)
    env["PYTHONPYCACHEPREFIX"] = str(cache / "pyc")
    return env


def start_sim(d: Path, gpu: int) -> subprocess.Popen:
    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    fh = (d / "run.log").open("w")
    return subprocess.Popen([str(EXE), "./params.txt"], cwd=str(d),
                            stdout=fh, stderr=subprocess.STDOUT,
                            env=env, start_new_session=True)


def free_gb(path: Path) -> float:
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize / 1e9


def cli_value(flag: str) -> str | None:
    """Value of `--flag value` or `--flag=value`, or None."""
    for i, a in enumerate(sys.argv):
        if a == flag and i + 1 < len(sys.argv):
            return sys.argv[i + 1]
        if a.startswith(flag + "="):
            return a.split("=", 1)[1]
    return None


def selected_arms() -> list[tuple[str, list[str]]]:
    """ARMS filtered by `--only name1,name2` (default: all)."""
    only = cli_value("--only")
    if not only:
        return list(ARMS)
    want = [s.strip() for s in only.split(",") if s.strip()]
    known = {name for name, _ in ARMS}
    unknown = [w for w in want if w not in known]
    if unknown:
        raise SystemExit(f"unknown arm(s): {unknown}; known: {sorted(known)}")
    return [(n, e) for n, e in ARMS if n in want]


def selected_gpus() -> list[int]:
    """GPUS filtered by `--gpus 4,6` (default: all). Lets a second wave use
    only the devices an earlier wave has freed."""
    spec = cli_value("--gpus")
    if not spec:
        return list(GPUS)
    return [int(s) for s in spec.split(",") if s.strip()]


def preflight(arms: list[tuple[str, list[str]]]) -> list[str]:
    """Blocking problems, as a list of human-readable strings."""
    problems = []
    if not BASE.is_file():
        problems.append(f"baseline params missing: {BASE}")
    if not EXE.is_file():
        problems.append(f"binary missing: {EXE}")
    else:
        try:
            blob = EXE.read_bytes()
            if b"rl_pump_kp_phantom" not in blob:
                problems.append(
                    f"{EXE.name} predates the per-sector gains (no "
                    f"rl_pump_kp_phantom symbol) -- rebuild before launching")
        except OSError as exc:
            problems.append(f"cannot read {EXE}: {exc}")
    if not VENV_PY.is_file():
        problems.append(f"venv python missing: {VENV_PY}")
    SCRATCH.mkdir(parents=True, exist_ok=True)
    need = max(MIN_FREE_GB * len(arms) / max(len(ARMS), 1), 20.0)
    avail = free_gb(SCRATCH)
    if avail < need:
        problems.append(
            f"only {avail:.0f} GB free on {SCRATCH}; need >= {need:.0f} GB "
            f"for {len(arms)} runs x {KEEP_LAST}+1 plotfiles")
    for name, _ in arms:
        if (OUT / name / "run.log").exists():
            problems.append(f"{name} already has a run.log -- would clobber; "
                            f"move it aside or pass --force")
    return problems


def main() -> int:
    global OUT
    out_spec = cli_value("--out")
    if out_spec:
        OUT = Path(out_spec) if "/" in out_spec else ROOT / "runs" / out_spec
    OUT.mkdir(parents=True, exist_ok=True)
    force = "--force" in sys.argv
    arms = selected_arms()
    gpus = selected_gpus()

    if "--dry-run" in sys.argv:
        for name, extra in arms:
            d = prepare(name, extra)
            print(f"{name}: {', '.join(extra)}\n   -> {d}")
        print("\n=== tail of generated params (first arm) ===")
        first = OUT / arms[0][0] / "params.txt"
        print("\n".join(first.read_text().splitlines()[-10:]))
        print("\n=== duplicate-key audit (must all be 1) ===")
        watched = sorted(STRIP_ALWAYS | {
            k for _, extra in ARMS for line in extra
            if (k := param_key(line)) is not None
        } | {"recipe_initial_data_file"})
        bad = 0
        for name, _ in arms:
            lines = (OUT / name / "params.txt").read_text().splitlines()
            counts = {k: sum(1 for ln in lines if param_key(ln) == k)
                      for k in watched}
            # A key only matters if this arm actually sets it.
            offenders = {k: c for k, c in counts.items() if c > 1}
            missing = [k for k in ("output_path", "amr.plot_file",
                                   "amr.check_file", "rl_pump_stop_time")
                       if counts[k] != 1]
            flag = ""
            if offenders or missing:
                bad += 1
                flag = f"  <-- DUPLICATE {offenders} MISSING {missing}"
            print(f"{name}: " + " ".join(f"{k}={counts[k]}" for k in watched
                                         if counts[k]) + flag)
        print("\n=== preflight ===")
        problems = preflight(arms)
        for p in problems:
            print(f"  BLOCK: {p}")
        if not problems:
            print("  all clear")
        print(f"  free on {SCRATCH}: {free_gb(SCRATCH):.0f} GB")
        print(f"  would use gpus: {gpus}")
        return 1 if bad else 0

    problems = [p for p in preflight(arms)
                if not (force and "run.log" in p)]
    if problems:
        for p in problems:
            log(f"[preflight-BLOCK] {p}")
        return 2

    log(f"confinement campaign A start: {len(arms)} arms "
        f"({', '.join(n for n, _ in arms)}) on gpus {gpus}, pump always on, "
        f"plot_interval={PLOT_INTERVAL}, "
        f"--delete --keep-last {KEEP_LAST}, {free_gb(SCRATCH):.0f} GB free")
    pending = list(arms)
    active: dict[int, dict] = {}

    while pending or active:
        for gpu in gpus:
            if gpu in active or not pending:
                continue
            name, extra = pending.pop(0)
            d = prepare(name, extra)
            cons = start_consumer(d)
            time.sleep(5)
            sim = start_sim(d, gpu)
            active[gpu] = {"name": name, "sim": sim, "consumer": cons, "dir": d}
            log(f"[start] {name} gpu={gpu} sim_pid={sim.pid} :: {'; '.join(extra)}")

        time.sleep(30)

        for gpu in list(active):
            job = active[gpu]

            # Abort a run that has already collapsed rather than spending the
            # GPU on it. The arm still drains, scores and keeps its data below
            # -- a killed arm is a measurement, not a lost run.
            if job["sim"].poll() is None:
                reason = collapse_verdict(job["dir"])
                if reason:
                    job["kill_reason"] = reason
                    log(f"[watchdog-KILL] {job['name']}: {reason}")
                    try:
                        os.killpg(os.getpgid(job["sim"].pid), signal.SIGTERM)
                        job["sim"].wait(timeout=120)
                    except subprocess.TimeoutExpired:
                        os.killpg(os.getpgid(job["sim"].pid), signal.SIGKILL)
                    except (ProcessLookupError, PermissionError) as exc:
                        log(f"[watchdog] {job['name']}: kill failed: {exc}")

            rc = job["sim"].poll()
            if rc is None:
                continue
            log(f"[sim-exit] {job['name']} rc={rc}; draining {DRAIN_SECONDS}s")
            time.sleep(DRAIN_SECONDS)
            cons = job["consumer"]
            if cons.poll() is None:
                os.killpg(os.getpgid(cons.pid), signal.SIGTERM)
                try:
                    cons.wait(timeout=180)
                except subprocess.TimeoutExpired:
                    os.killpg(os.getpgid(cons.pid), signal.SIGKILL)

            # governor_min plus the time it first crosses 0.5: an arm that
            # chokes before t=27 bought its confinement with constraint
            # violation and is disqualified.
            gov_min, gov_cross = "?", "never"
            cn = job["dir"] / "data/constraint_norms.dat"
            if cn.exists():
                try:
                    vals = []
                    for ln in cn.read_text().splitlines():
                        f = ln.split()
                        if len(f) > 9:
                            vals.append((float(f[0]), float(f[9])))
                    if vals:
                        gov_min = f"{min(v for _, v in vals):.4f}"
                        below = [t for t, v in vals if v < 0.5]
                        if below:
                            gov_cross = f"t={below[0]:.2f}"
                except (ValueError, IndexError):
                    pass

            log(f"[score] {job['name']}: 4D evolving geodesic")
            try:
                sc_rc = subprocess.run(
                    [str(VENV_PY), "-u", str(SCORE_EVOLVING), str(job["dir"])],
                    env=consumer_env(), capture_output=True, text=True,
                    timeout=3600,
                )
                for line in (sc_rc.stdout or "").splitlines():
                    if line.strip():
                        log(f"    {line.rstrip()}")
                if sc_rc.returncode != 0:
                    log(f"[score-FAIL] {job['name']} rc={sc_rc.returncode} "
                        f"{(sc_rc.stderr or '').strip()[:300]}")
            except subprocess.TimeoutExpired:
                log(f"[score-FAIL] {job['name']}: timed out after 3600s")
            except Exception as exc:  # noqa: BLE001
                log(f"[score-FAIL] {job['name']}: {exc}")

            # Scratch is transient: anything not extracted before the purge is
            # gone for good. Purge ONLY when every resident plotfile is
            # recorded as successfully extracted.
            sc = scratch_dir(job["name"])
            leftover = list(sc.glob("RadialRecipePlt*"))
            state_f = job["dir"] / "small_data/consume_state.json"
            extracted = set()
            try:
                extracted = set(json.loads(state_f.read_text()))
            except (OSError, ValueError):
                pass
            unextracted = [p.name for p in leftover if p.name not in extracted]
            if unextracted:
                log(f"[KEEP-SCRATCH] {job['name']}: {len(unextracted)} "
                    f"plotfile(s) NOT extracted, keeping {sc} -- {unextracted}")
            else:
                shutil.rmtree(sc, ignore_errors=True)
            killed = job.get("kill_reason")
            log(f"[done] {job['name']} rc={rc} governor_min={gov_min} "
                f"governor_cross={gov_cross} leftover_plt={len(leftover)} "
                f"free={free_gb(SCRATCH):.0f}GB"
                + (f" watchdog={killed}" if killed else ""))
            del active[gpu]

    log("confinement campaign A complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
