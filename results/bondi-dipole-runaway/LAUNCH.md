# Launch reference

Exact commands and configuration for every cell in the pack. Paths are
repo-relative; run from the repository root.

---

## 1. Prerequisites

```bash
# Constraint solver (complex scalar matter, per-lump signs)
bash grteclyn-wrapper/scripts/wormhole/build/build_grtresna_bosonstar.sh
# -> GRTresna/Examples/BosonStarBH/Main_BosonStarBH3d.*.ex
```

Needs `CHOMBO_HOME` and the solver toolchain on `PATH` (see
`grteclyn-wrapper/README.md`). The GRTeclyn `RadialRecipe` executable is built
per that README as well.

**Single rank everywhere.** MPI is unusable on this node (`mpirun` segfaults),
so both the solve and the evolution run with one rank; the launchers hard-set
`GRTRESNA_RANKS=1` and must not be raised without re-testing.

## 2. The six cells

Each launcher registers a stop handle (`launcher.pid`) and writes to
`runs/bondi/<cell>/`. All six share the same physics configuration — only the
sector flags, the per-lump frequency, and the stop time differ.

```bash
# --- calibration singles (stop 40) ------------------------------------------
BONDI_GPU=1 BONDI_RUNS_DIR="$PWD/runs/bondi/single_p" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_single_selfgrav.sh

BONDI_EXOTIC=1 BONDI_GPU=2 BONDI_RUNS_DIR="$PWD/runs/bondi/single_m" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_single_selfgrav.sh

# --- the runaway cell (stop 60) ---------------------------------------------
BONDI_S0=0 BONDI_S1=1 BONDI_GPU=3 BONDI_RUNS_DIR="$PWD/runs/bondi/pair_pm" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh

# --- controls (stop 60) -----------------------------------------------------
BONDI_S0=0 BONDI_S1=0 BONDI_GPU=1 BONDI_RUNS_DIR="$PWD/runs/bondi/pair_pp" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh

BONDI_S0=1 BONDI_S1=1 BONDI_GPU=2 BONDI_RUNS_DIR="$PWD/runs/bondi/pair_mm" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh

# --- equal-|ADM| variant: phantom at its own frequency (stop 60) ------------
BONDI_S0=0 BONDI_S1=1 BONDI_S1_OMEGA=0.56598 BONDI_GPU=0 \
  BONDI_RUNS_DIR="$PWD/runs/bondi/pair_pm_eqm" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh
```

Detached (the campaign ran this way, one cell per GPU):

```bash
setsid nohup /usr/bin/env BONDI_S0=0 BONDI_S1=1 BONDI_GPU=3 \
  BONDI_RUNS_DIR="$PWD/runs/bondi/pair_pm" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh \
  > runs/bondi/pair_pm.launch.log 2>&1 < /dev/null & disown
```

> **`/usr/bin/env`, spelled out — not a style preference.** On the GPU build node
> two other users' bin directories precede `/usr/bin` on `PATH` (check with
> `type -a env`), and each contains an executable named `env` that is really a
> PATH-setup snippet meant to be *sourced*. Run as `env VAR=x cmd` it prepends
> its own directory to `PATH`
> and **exits 0 without ever executing `cmd`**. A bare `env` here therefore
> launches nothing at all, leaving an empty `.launch.log` and a success exit
> code — indistinguishable at a glance from a launch that is merely slow to
> write output. (Verified 2026-08-17; only `env` is shadowed, no other core tool.)

Verify a detached launch with `pgrep -f "bondi_sg_pair_p[m]"` before concluding
anything from the log — the log file can lag the detach, and per the note above
a silent no-op looks exactly the same.

**Stopping** is only ever done with the sanctioned tool, which kills the
orchestrator first and then sweeps workers:

```bash
bash grteclyn-wrapper/scripts/campaigns/stop_campaign.sh runs/bondi/<cell>
```

Killing a worker directly just advances the orchestrator to the next cell.

### Launcher knobs

| variable | `run_single_selfgrav.sh` | `run_pair_selfgrav.sh` |
|---|---|---|
| `BONDI_GPU` | default 1 | default 3 |
| `BONDI_STOP_TIME` | default 40 | default 60 |
| `BONDI_SEP` | 8 (R₀ = 4) | 8 (R₀ = 4) |
| `BONDI_RUNS_DIR` | output parent | output parent |
| `BONDI_EXOTIC` | `1` ⇒ lone phantom (`bondi_sg_single_m`) | — |
| `BONDI_S0` / `BONDI_S1` | — | sector per lump: 0 canonical, 1 phantom ⇒ suffix `pm`/`pp`/`mm` |
| `BONDI_S1_OMEGA` | — | per-lump star frequency for lump 1; appends `_eqm` to the run name |

### Rerun knobs (added 2026-08-17 for Debug.md §3)

Every one defaults to the published campaign value, so an un-set environment
reproduces the original cell. Both launchers accept all of them.

| variable | default | what it is for |
|---|---|---|
| `BONDI_NFULL` | 128 | item A, convergence series (192, 256) |
| `BONDI_MAXLEVEL` | 1 | item A: **set 0** for convergence cells — see the warning below |
| `BONDI_LFULL` | 64 | item C, enlarged domain (128, 256) |
| `BONDI_RADII` | `8 16` | item C, wave-zone extraction shells |
| `BONDI_DT_MULT` | 0.02 | item I, CFL experiment — **do not raise, see below** |
| `BONDI_SCRUTINY` | 0 | item B, momentum-balance stream (`sector_dynamics.dat`) |
| `BONDI_PSI4_HIGHER_L` | 0 | item C, l≥3 Psi4 stream (`psi4_mode_higher_l.dat`) |
| `BONDI_PSI4_ELLS` | `3 4` | which multipoles that stream carries |

**Two measured warnings.**

* `BONDI_DT_MULT=0.2` **fails** — the star disperses (rms 5.05 → 19.2 by t = 40
  against a ±10 % gate) and the tagger chases the resulting noise until the
  refined region is 380× larger. Keep 0.02. Full numbers in Debug.md item I.
* Convergence cells need `BONDI_MAXLEVEL=0`, for memory as much as for physics.
  A refined N = 128 cell measured **41 GB** of GPU memory against 8.8 GB
  unigrid; carried to N = 256 that would OOM an 80 GB H100 outright. Unigrid is
  also what Richardson extrapolation requires, and the t ≤ 30 science window was
  effectively unigrid anyway.

Memory to budget, measured on this node: **8.8 GB per unigrid N = 128 cell**,
scaling as N³ — so roughly 30 GB at N = 192 and 55–60 GB at N = 256. One
N = 256 cell per GPU, and do not co-schedule anything with it.

Why the stop times: singles are calibration (the pass gate is rms flat ±10 %
to t = 40 with min χ plateau above 0.3); mixed pairs are capped at **60**
because the massive-field radiation bath cannot exit the massless-wave
boundaries and begins to contaminate the drift measurement beyond that.

## 3. What the launchers set

Both scripts replay a stored evaluation for the numerics and override the
physics. The full expanded inputs are packed per cell as
`data/<cell>/grtresna_params.txt` and `data/<cell>/evolution_params.txt`.

**Grid / evolution**

| setting | value |
|---|---|
| domain, resolution | `L_full = 64`, `N_full = 128`, centre `32 32 32` |
| AMR | `max_level = 1`, `regrid_threshold = 0.02`, `regrid_interval = 16` |
| time step | `dt_multiplier = 0.02` |
| boundaries | Sommerfeld (`hi_boundary = lo_boundary = 1 1 1`), non-periodic |
| plot cadence | `plot_interval = 40` (152 frames over t = 0–60) |
| extraction | radii 8 and 16 |

**Constraint solve**

| setting | value |
|---|---|
| solver | GRTresna CTTKHybrid, `BosonStarBH`, `--grtresna-max-level 3`, `--grtresna-domain-l 128` |
| iterations / tolerances | 50 max, `NL_exit_tolerance = 0.1 %`, `NL_stall_tolerance = 0.002` |
| ranks | 1 |

**Matter overrides (identical across cells)**

```
grtresna_scalar_lambda = 10240        grtresna_scalar_mu = 21845333
grtresna_bs_omega      = 0.55         grtresna_bs_selfgrav = 1
trajectory_well_width  = 1.2          trajectory_lump{k}_well_depth = 0.15
rl_pump_stop_time      = 0            grtresna_boost_lumps = 0
trajectory_A_breath = trajectory_z_amp = trajectory_omega_z = 0
trajectory_lump{k}_{R0=4, tilt_theta=0, tilt_phi=0, v_rad=0, omega_rot=0}
trajectory_lump0_phase0 = 0           trajectory_lump1_phase0 = pi
trajectory_lump{k}_exotic = <sector>  trajectory_lump1_bs_omega = <eqm only>
```

Lumps start **at rest**: no orbital or radial velocity, no boost in the initial
data, no breathing, no tilt, no rotation, and the trajectory pump is off for the
entire evolution (`rl_pump_stop_time = 0`). Any drift is therefore
gravitational, not driven.

**Diagnostics** are switched on by the launchers:

```bash
export GRTECLYN_SECTOR_BARYCENTERS=1   # small_data/sector_barycenters.dat
export GRTECLYN_PSI4=1                 # small_data/psi4_*.dat
```

## 4. Verifying a launch before trusting it

Check the t = 0 row of `small_data/sector_barycenters.dat` against the star
predictions — this catches a mis-seeded run immediately:

| seed | total | rms |
|---|---|---|
| canonical star, ω = 0.550 | 15.92 | 5.05 |
| phantom star, ω = 0.550 | 20.99 | 5.43 |
| phantom star, ω = 0.56598 | 17.05 | 5.16 |

Pair cells with one lump per sector show both fingerprints in the same row;
same-sector pairs show twice the single value. Also confirm the final row of
`grtresna/Ham_and_Mom_errors.txt` is inside the 0.1 % gate, and that
`grtresna/params.txt` carries `lump0_amp ≈ 0.019695` (the solved φ_c, not a
clamped amplitude) — the failure mode described in [`DEBUGGING.md`](DEBUGGING.md) §1.

## 5. Post-processing

```bash
# Movies from the retained PNG frames (19 fields, or a subset with --only)
bash grteclyn-wrapper/scripts/plot/make_movies.sh runs/bondi/pair_pm/bondi_sg_pair_pm \
     --framerate 10 --only scalar_activity_proj_z chi_minus_1_z rho_req_z

# Rebuild this results pack (light artefacts only, machine paths scrubbed)
bash research/bondi_dipole/pack_results.sh

# Derived tables
python3 results/bondi-dipole-runaway/analysis/make_tables.py
python3 results/bondi-dipole-runaway/analysis/newtonian_reference.py
PYTHONPATH=grteclyn-wrapper/src grteclyn-wrapper/.venv/bin/python \
  results/bondi-dipole-runaway/analysis/star_family_scan.py
```

## 6. Cost

Each pair cell is one constraint solve (a few minutes, single rank) plus a
GPU evolution to t = 60 at N = 128 — roughly 2–3 hours per cell on one device;
singles to t = 40 are proportionally shorter. The four-cell control set and the
equal-mass variant were run concurrently, one per GPU.
