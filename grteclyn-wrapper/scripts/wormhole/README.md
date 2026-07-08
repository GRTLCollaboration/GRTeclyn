# Rotating Wormhole Scripts

Pipeline for constraint-clean rotating-wormhole evolutions: build the solver,
generate initial data, run the evolution with disk-safe frame streaming, and
post-process.

```
 build/     id/              run/                 postrun/
 GRTresna   solve kappa      evolve + frames      analysis / archive
 solver     family + Teo     (sidecar deletes
            analytic          plotfiles on the fly)
```

## Workflow

```bash
# 0. Build the GRTresna constraint solver (once)
bash scripts/wormhole/build/build_grtresna_bosonstar.sh

# 1. Solve the kappa initial-data family (constraint-clean .gridinit files)
RES_N=128 bash scripts/wormhole/id/solve_kappa_family.sh 1.0,0.7,0.5 4

# 2. Evolve a single case (params generated from template, frames streamed)
bash scripts/wormhole/run/wormhole_case.sh \
    --kappa 1.0 --dx 0.5 --max-level 3 --gpu 0

# 2b. Larger box for boundary convergence study
EVO_L=96 RES_N=192 bash scripts/wormhole/id/solve_kappa_family.sh 1.0 4
bash scripts/wormhole/run/wormhole_case.sh \
    --kappa 1.0 --dx 0.5 --box-size 96 --gpu 0
```

## Directory layout

### `build/`

| Script | Purpose |
|--------|---------|
| `build_grtresna_bosonstar.sh` | Rebuild the GRTresna `BosonStarBH` MPI solver (`CTTKHybrid<ComplexScalarField>`). |

### `id/` — initial data

| Script | Purpose |
|--------|---------|
| `solve_kappa_family.sh` | Shell wrapper for the kappa-family solver. |
| `solve_kappa_family.py` | Solve Ham + Mom constraints for each kappa, export `.gridinit` at evolution dx. Env: `EVO_L` (box size, default 64), `RES_N` (grid cells, default 64). |
| `make_teo_wormhole_gridinit.py` | Analytic Teo rotating wormhole `.gridinit` (weak-spin baseline, not constraint-solved). |

### `run/` — evolution

| Script | Purpose |
|--------|---------|
| `wormhole_case.sh` | CLI evolution driver (shell wrapper). |
| `wormhole_case.py` | Generate params from template + CLI flags, launch GRTeclyn with frame/Psi4 consumer sidecar (L6 disk discipline). Key flags: `--kappa --dx --box-size --max-level --stop-time --gpu`. |
| `run_rotating_wormhole.sh` | Generic launcher: run any `params_*.txt` with the consumer sidecar. |

### `postrun/` — analysis and utilities

| Script | Purpose |
|--------|---------|
| `scan_rotating_wormhole_support.py` | Sweep runtime phantom support strength across cases. |
| `move_files.sh` | Archive run outputs (frames, diagnostics) to `output/SimResults/`. |
| `rollback` | Roll back a run directory to a checkpoint index N (delete later data). |

## Key conventions

- **L2 (critical):** never evolve finer than the `.gridinit`'s native dx.
  `--dx` selects both the GRTresna solve resolution and the evolution level-0 dx.
- **L6 (disk discipline):** the consumer sidecar streams frames + Psi4 small_data
  and deletes plotfiles during the run. A finished run is ~15 MB, not tens of GB.
- **Stop-time auto-scales:** when `--stop-time` is not set, defaults to
  `r_outer + 6` (light-crossing to the outermost extraction shell + buffer).
- **Extraction:** 2 radii (inner=12, outer=L/2-8) to keep the consumer lightweight.
- **Tags:** directory names encode `(omega, m, kappa, dx)` and append `_L<N>`
  only when the box size differs from the default 64 (backward compatible).
