#!/usr/bin/env python3
"""Clean pump dose-response ladder, controller_reservoir_mode = 0 throughout.

Why mode 0: the reservoir ledger diverges (pressureless dust has no restoring
force), and in mode >= 1 it also fed the pump safety governor, which throttled
the pump to zero at t~10 in every run of the first campaign. Mode 0 removes
both problems. Constraint defence comes from the Duhamel bound
(pump_constraint_budget.py) plus the pump-free comparison rung.

All five runs clone the same RM baseline params and share the cached gridinit,
so only rl_pump_stop_time varies. tp4 doubles as a bit-identity regression
check against the earlier runs/pump/always_on_pump/hq146_m0_tp4_t30 -- the governor
fix must be a no-op for mode 0.

Requires binary d6a0c350 or later.

REFERENCE IMPLEMENTATION of the mandatory node-local plotfile scratch pattern
(see scripts/campaigns/README.md, "Plotfile scratch is node-local"). The parts
worth copying into any new campaign launcher:

  * amr.plot_file / amr.check_file -> SCRATCH (local NVMe); output_path stays
    on NFS, so only .dat and small_data/ cross the network.
  * the baseline clone STRIPS amr.plot_file / amr.check_file, which are
    absolute paths into the source run dir and would otherwise let a --delete
    consumer prune the baseline's own plotfiles.
  * the consumer runs `.venv/bin/python -m ...` rather than `uv run`, with
    every cache/temp env var pinned inside SCRATCH.
  * scratch is purged only when every resident plotfile appears in
    consume_state.json; otherwise it is kept and logged [KEEP-SCRATCH].
  * DRAIN_SECONDS = 600 after sim exit before the consumer is torn down.

ROOT/BASE/GPUS below are specific to the local cluster deployment -- adjust them,
keep the structure.
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

# grteclyn-wrapper/scripts/campaigns/rl/ → repo root
ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / "runs/grtresna_promote/bcma_rm_L128_N256_t30_hq_eval000146/params.txt"
OUT = ROOT / "runs/pump/pump_ladder_m0"
EXE = ROOT / "Examples/RadialRecipe/main3d.gnu.CUDA.ex"
WRAPPER = ROOT / "grteclyn-wrapper"
# Call the project venv python DIRECTLY rather than via `uv run`.
# uv writes to ~/.cache/uv; going straight to the venv keeps every
# write inside the project, /tmp, or the NFS run dir. Nothing this
# queue does touches ~/.local (which is admin-locked anyway).
VENV_PY = WRAPPER / ".venv/bin/python"
SCORE_EVOLVING = WRAPPER / "scripts/campaigns/rl/score_evolving_geodesic.py"

# Plotfiles are 3.2 GB each and are written AND read back every ~288 s per run.
# On NFS that is ~22 MB/s per run plus heavy metadata churn, which is what
# caps concurrency. amr.plot_file is independent of output_path, so the big
# transient data goes to node-local NVMe (the overlay mount) and only the
# small .dat / small_data results ever touch NFS.
SCRATCH = Path("/tmp/grteclyn_scratch")

MAX_CONCURRENT = 6
GPUS = [0, 1, 2, 3, 4, 5]
PLOT_INTERVAL = 144
DRAIN_SECONDS = 600

# (name, rl_pump_stop_time)   -1.0 == pump never stops
MATRIX = [
    ("lad_m0_tp30", -1.0),  # always on -- headline + lapse-collapse test
    ("lad_m0_tp4", 4.0),    # regression check vs the earlier m0_tp4
    ("lad_m0_tp16", 16.0),  # last rung that stayed healthy at fast tier
    ("lad_m0_tp24", 24.0),  # brackets the turnover between 16 and 30
    ("lad_m0_tp8", 8.0),
    ("lad_m0_tp0", 0.0),    # pump-free rung
]


def log(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with (OUT / "queue.log").open("a") as fh:
        fh.write(line + "\n")


def scratch_dir(name: str) -> Path:
    return SCRATCH / name


def prepare(name: str, tstop: float) -> Path:
    d = OUT / name
    (d / "data").mkdir(parents=True, exist_ok=True)
    scratch = scratch_dir(name)
    scratch.mkdir(parents=True, exist_ok=True)
    keep = []
    for line in BASE.read_text().splitlines():
        s = line.strip()
        # amr.plot_file / amr.check_file are absolute paths into the RM
        # baseline dir and MUST be redirected.
        if s.startswith(("output_path", "plot_interval", "rl_pump_stop_time",
                         "controller_reservoir_mode", "amr.plot_file",
                         "amr.check_file")):
            continue
        keep.append(line)
    keep += [
        "",
        "# --- clean mode-0 pump dose-response ladder ---",
        f'output_path = "{d}"',
        f'amr.check_file = "{scratch}/RadialRecipeChk"',
        f'amr.plot_file  = "{scratch}/RadialRecipePlt"',
        f"plot_interval       = {PLOT_INTERVAL}",
        f"rl_pump_stop_time = {tstop}",
        "controller_reservoir_mode = 0",
    ]
    (d / "params.txt").write_text("\n".join(keep) + "\n")
    return d


def start_consumer(d: Path) -> subprocess.Popen:
    cmd = [
        str(VENV_PY), "-m",
        "grteclyn_wrapper.visualisation.process_wave.consume_plotfiles",
        "--data", str(scratch_dir(d.name)), "--out", str(d / "small_data"),
        "--center", "64", "64", "64",
        "--radii", "8", "12", "24", "--n-points", "32",
        "--no-psi4", "--ftl-timeseries", "--ftl-l", "8",
        "--evolving-geodesic",
        "--confinement-timeseries", "--confinement-well-width", "1.667",
        "--watch", "--delete", "--keep-last", "3",
        "--stable-seconds", "30", "--poll-seconds", "20", "-j", "1", "--verbose",
    ]
    env = consumer_env()
    fh = (d / "consumer.log").open("w")
    return subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT,
                            env=env, start_new_session=True)


def consumer_env() -> dict:
    """Environment shared by the consumer AND the post-run scoring pass.

    They must agree: the scoring pass reads the cache the consumer wrote, so a
    mismatch in GRTECLYN_METRIC_STACK_N_SPACE would have it validate against
    the wrong resolution.
    """
    env = dict(os.environ)
    env["GEODESIC_EMIT_MIN_TIME"] = "0"
    env.pop("RL_PUMP_STOP_TIME", None)
    # metric_stack resolution MUST match the run's finest AMR level.  The cache
    # is a UNIFORM resample: at n_space=33 its dx is 0.5, the COARSEST level,
    # so every feature living on the 3 refined levels (finest dx=0.0625) is
    # smoothed away with no error and no artifact.  On the 2026-07-28 ladder
    # this erased the pump-free run's collapsing well (true min_chi 5.7e-4,
    # cached 5.6e-2 -- 99x too shallow); its rays never paid the Shapiro delay
    # and it reported the LARGEST apparent shortcut precisely because it was
    # the worst-resolved run.  257 = 2*half_width/finest_dx + 1.
    # See metric_stack_cache.required_n_space / cache_fidelity.
    env["GRTECLYN_EVOLVING_GEODESIC_MODE"] = "hq"
    env["GRTECLYN_METRIC_STACK_N_SPACE"] = "257"
    # Keep every cache/temp write inside our own scratch, never ~/.local
    # or ~/.cache.
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


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    if "--dry-run" in sys.argv:
        for name, tstop in MATRIX:
            d = prepare(name, tstop)
            print(f"{name} (tpump={tstop}) -> {d}")
        d = OUT / MATRIX[0][0]
        print("\n=== tail of generated params ===")
        print("\n".join((d / "params.txt").read_text().splitlines()[-9:]))
        for key in ("output_path", "amr.plot_file", "amr.check_file",
                    "recipe_initial_data_file", "rl_pump_stop_time",
                    "controller_reservoir_mode", "plot_interval"):
            hits = [ln for ln in (d / "params.txt").read_text().splitlines()
                    if ln.strip().startswith(key)]
            print(f"{key}: {len(hits)} occurrence(s)")
        return 0

    log(f"ladder start: {len(MATRIX)} runs, mode 0, max {MAX_CONCURRENT} "
        f"concurrent, plot_interval={PLOT_INTERVAL}, --delete --keep-last 3")
    pending = list(MATRIX)
    active: dict[int, dict] = {}

    while pending or active:
        for gpu in GPUS:
            if gpu in active or not pending:
                continue
            name, tstop = pending.pop(0)
            d = prepare(name, tstop)
            cons = start_consumer(d)
            time.sleep(5)
            sim = start_sim(d, gpu)
            active[gpu] = {"name": name, "sim": sim, "consumer": cons, "dir": d}
            log(f"[start] {name} gpu={gpu} tpump={tstop} sim_pid={sim.pid}")

        time.sleep(30)

        for gpu in list(active):
            job = active[gpu]
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
            gov_min = "?"
            cn = job["dir"] / "data/constraint_norms.dat"
            if cn.exists():
                try:
                    vals = [float(l.split()[9]) for l in
                            cn.read_text().splitlines() if len(l.split()) > 9]
                    gov_min = f"{min(vals):.4f}" if vals else "?"
                except (ValueError, IndexError):
                    pass
            # The 4D evolving trace is NOT run by the consumer.  --evolving-geodesic
            # only builds the metric_stack cache; the trace itself lives in the
            # metrics aggregation layer, which only the QD/CMA-ES evaluation path
            # calls.  Without this step f_geo_evol / f_geo_evol_ok (cols 13/14 of
            # ftl_timeseries.dat) keep the hardcoded "0.0  0" placeholders that
            # extraction/ftl.py writes on every row -- they look computed-and-zero
            # rather than never-computed, which is exactly how that went unnoticed
            # for two campaigns.  Runs BEFORE the scratch purge below so the
            # plotfile fallback is still available if the cache is unusable.
            log(f"[score] {job['name']}: 4D evolving geodesic")
            try:
                sc_rc = subprocess.run(
                    [str(VENV_PY), "-u", str(SCORE_EVOLVING), str(job["dir"])],
                    env=consumer_env(), capture_output=True, text=True, timeout=3600,
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

            # Scratch is transient: anything not extracted before we purge is
            # gone for good (on NFS it would have survived). So purge ONLY when
            # every resident plotfile is recorded as successfully extracted.
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
                log(f"[KEEP-SCRATCH] {job['name']}: {len(unextracted)} plotfile(s) "
                    f"NOT extracted, keeping {sc} -- {unextracted}")
            else:
                shutil.rmtree(sc, ignore_errors=True)
            log(f"[done] {job['name']} rc={rc} governor_min={gov_min} "
                f"leftover_plt={len(leftover)}")
            del active[gpu]

    log("ladder complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
