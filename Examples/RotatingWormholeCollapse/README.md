# RotatingWormholeCollapse

Evolution of a **rotating traversable wormhole** (complex phantom scalar,
Ellis–Bronnikov / Kleihaus–Kunz class) in GRTeclyn. Rotation supplies an inherent
quadrupole, so the ℓ=2 gravitational-wave signal is *natural* (no artificial
perturbation). The rotating ID is produced by a **constraint-clean GRTresna
solve**, which removes the O(ω) momentum-constraint defect that made naive
analytic initial data diverge.

Full program, status, and physics rationale:
[`grteclyn-wrapper/docs/RotatingWormholePlan.md`](../../grteclyn-wrapper/docs/RotatingWormholePlan.md).

---

## Pipeline (3 stages)

```
 (1) GRTresna solve            (2) flatten to .gridinit         (3) GRTeclyn evolve + frames
 exotic complex winding   ->   uniform grid at evolution   ->   complex_scalar model, CCZ4,
 scalar on a singular-ψ        level-0 dx (L2-safe)             Weyl4 extraction, diagnostics
 throat (BosonStarBH)          convert_chombo_to_gridinit       main3d.gnu.MPI.CUDA.ex
```

Both stages are driven by **CLI scripts** in
`grteclyn-wrapper/scripts/wormhole/{id,run}/` — params are generated from
templates, so we do **not** commit a `params_*.txt` per case (same philosophy as
the QD campaign). See `grteclyn-wrapper/scripts/wormhole/README.md` for the full
directory layout.

- **Matter model** `complex_scalar`: Φ = f(r,θ) e^{i(mφ_az)} stored as
  (φ, Π, φ2, Π2). |Φ|² is axisymmetric (no four-lobe), J_z carried by the phase
  winding. Phantom (−ρ) supports the throat.
- **κ (amplitude-reduction) trigger**: solve the constraints for field amplitude
  `f → κ·f`. κ=1 is the equilibrium; κ<1 is an *exact* Einstein solution out of
  equilibrium that collapses — with **zero** initial constraint defect.
- **Lesson L2 (critical)**: never evolve finer than the ID file's native dx.
  `--dx` selects both the GRTresna solve resolution and the evolution level-0 dx:
  `dx=1.0 → N=64`, `dx=0.5 → N=128`.

---

## Main commands

All commands run from `grteclyn-wrapper/` (they source `scripts/lib/env.sh`).

### 0. Build the binaries
```bash
# GRTeclyn evolution (GPU)
make -C Examples/RotatingWormholeCollapse -j8 USE_CUDA=TRUE USE_MPI=TRUE COMP=gnu CUDA_ARCH=90
# GRTresna initial-data solver (CPU/MPI)
bash grteclyn-wrapper/scripts/wormhole/build/build_grtresna_bosonstar.sh
```

### 1. Solve the κ initial-data family (produces `.gridinit`s)
```bash
# dx=0.5 (N=128) family, kappas 1.0/0.7/0.5, 4 MPI ranks:
RES_N=128 bash grteclyn-wrapper/scripts/wormhole/id/solve_kappa_family.sh 1.0,0.7,0.5 4
# dx=1.0 (N=64) quick family:
RES_N=64  bash grteclyn-wrapper/scripts/wormhole/id/solve_kappa_family.sh 1.0,0.9,0.7,0.5 2
```
Writes `runs/rotating_wormhole_id/rotwh_omega_p0p05_m1_kappa_<κ>_dx<dx>/initial_data.gridinit`
(the heavy solve HDF5s are auto-pruned; set `KEEP_SOLVE_SCRATCH=1` to keep them).
Reports Ham/Mom % per κ.

### 2. Evolve one case (params generated from template, frames rendered)
```bash
# high-res equilibrium arm:
bash grteclyn-wrapper/scripts/wormhole/run/wormhole_case.sh \
    --kappa 1.0 --dx 0.5 --max-level 3 --gpu 0
# collapse arms:
bash grteclyn-wrapper/scripts/wormhole/run/wormhole_case.sh --kappa 0.7 --dx 0.5 --gpu 1
bash grteclyn-wrapper/scripts/wormhole/run/wormhole_case.sh --kappa 0.5 --dx 0.5 --gpu 2

# boundary convergence study (larger box, boundary farther from throat):
EVO_L=96 RES_N=192 bash grteclyn-wrapper/scripts/wormhole/id/solve_kappa_family.sh 1.0 4
bash grteclyn-wrapper/scripts/wormhole/run/wormhole_case.sh \
    --kappa 1.0 --dx 0.5 --box-size 96 --gpu 0
```
Key flags: `--kappa --dx {0.5,1.0} --box-size --omega --m --max-level
--stop-time --plot-interval --gpu --no-frames --keep-plotfiles --dry-run`.
Stop-time auto-scales with extraction radius when not set explicitly
(`r_outer + 6`). Extraction uses 2 radii (inner=12, outer=L/2-8).
Output: `runs/rotating_wormhole/evo_omega_p0p05_m1_kappa_<κ>_dx<dx>_ml<L>/output/`
with `data/constraint_norms.dat`, `data/collapse_diagnostics.dat`, and
`frames/<field>_z/frames/*.png`. Plotfiles are pruned after frame extraction
(disk discipline L6); pass `--keep-plotfiles` to retain them.

### Diagnostics to check
- `data/constraint_norms.dat` — columns `t  Ham_L2  Mom_L2` (should stay bounded).
- `data/collapse_diagnostics.dat` — `t min_lapse min_chi ... J_z ...` (J_z conserved
  = spin held; min_chi→0 / horizon = collapse).

---

## Legacy analytic params (kept for reference / regression)

These pre-date the CLI and use the real-scalar (`exotic_scalar`) or analytic
per-level ID paths; the CLI above supersedes them for production:

- `params_rotating_complex_equilibrium.txt` — analytic complex-scalar equilibrium.
- `params_rotating_complex_collapse.txt` — analytic collapse (S_support trigger).
- `params_rotating_grtresna_{smoke,exotic,hires,a0_exotic,no_matter}.txt` —
  early GRTresna-ID + real `exotic_scalar` evolution tests.

Run any of them directly:
```bash
cd Examples/RotatingWormholeCollapse && CUDA_VISIBLE_DEVICES=0 \
  ./main3d.gnu.MPI.CUDA.ex params_rotating_grtresna_exotic.txt
```

---

## Status (2026-07-08)

Phase B complete: the GRTresna-solved rotating equilibrium **holds** where naive
analytic ID diverged.

| κ (dx=1.0, unigrid, t=8) | min_chi | min_lapse | max Ham | J_z | horizon |
|--------------------------|---------|-----------|---------|-----|---------|
| 1.0 (equilibrium) | 0.96 | 0.996 | 5.7e-3 | conserved | none |
| 0.7 | 0.86 | 0.963 | 4.6e-3 | conserved | none |
| 0.5 | 0.82 | 0.914 | 4.2e-3 | conserved | none |

Monotonic collapse trend with decreasing κ; constraints bounded; spin conserved.
Phase C (high-res N=128/dx=0.5, max_level=3, (ω,m,κ) grid) is in progress, then
Phase D (Ψ₄ multipole GW extraction + Kerr-QNM ringdown fit).
