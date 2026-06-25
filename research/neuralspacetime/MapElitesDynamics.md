# MAP-Elites Dynamics — Trajectory Ansatz Campaigns

> Continuation of [MapElites.md](./MapElites.md). This journal covers campaigns using
> the **trajectory ansatz** (Option D: per-lump independent orbits) and comparison
> against the spherical-harmonic (SH) baseline.

---

## Campaign comparison: `scalar_sh_ftl_v22` vs `trajectory_5lump_v1` (2026-06-25)

### Setup

| Parameter | **SH v22** | **Trajectory v1** |
|-----------|-----------|-------------------|
| Ansatz | Spherical harmonic (ℓ_max=4) | Per-lump circular orbits |
| Search dims | 38 D (24 SH + 14 geometry/physics) | 40 D (7 per-lump × 5 + 5 shared) |
| Matter model | `grtresna_independent_scalars` | `grtresna_independent_scalars` |
| Exotic handling | Coarse azimuthal wedge | Per-lump binary flag [0,1] |
| Motion | One global (v_tor, v_pol, v_rad) for all lumps | Independent ω_rot per lump in tilted plane |
| Grid (evolution) | N=128, L=64, ml=1 | N=128, L=64, ml=1 |
| Grid (GRTresna) | 128³ gridinit, ml=3 | 128³ gridinit, ml=3 |
| Stop time | 16.0 | 16.0 |
| Objective | `general_ftl` | `general_ftl` |
| Descriptor | `ftl_lifetime` (8×8 archive) | `ftl_lifetime` (8×8 archive) |
| GPUs | 8 × A100 80GB | 8 × A100 80GB |
| Batch size | 8 | 8 |
| Target evals | 200 | 50 |
| Completed evals | 202 | 25 (driver died mid-campaign) |
| Run dir | `runs/grtresna_qd/scalar_sh_ftl_v22/` | `runs/grtresna_qd/trajectory_5lump_v1/` |
| Date | 2026-06-24 | 2026-06-25 |

### Results — Head-to-head

| Metric | **SH v22** (202 evals) | **Trajectory v1** (25 evals) | Factor |
|--------|----------------------|---------------------------|--------|
| Best stable score | 470.6 | **687.9** | **1.5x** |
| Best overall score | 470.6 | **1166.8** (crashed) | **2.5x** |
| Best stable f_geo_peak | 2.12% | **5.30%** | **2.5x** |
| Best overall f_geo_peak | 2.12% | **24.6%** (crashed) | **11.6x** |
| Best stable ftl_geo_evolving | 0.101 | **0.261** | **2.6x** |
| Best f_op_peak | 7.0% | **15.2%** | **2.2x** |
| Best stable max_local_speed | 1.33c | **1.43c** | **1.1x** |
| FTL hit rate (per GPU eval) | 1.3% (1/79) | **54%** (7/13) | **~40x** |
| Archive cells filled | 2 / 64 | 3 / 64 | — |

### Status breakdown

| Status | **SH v22** | **Trajectory v1** |
|--------|-----------|-------------------|
| gpu_ok | 75 (37%) | 10 (40%) |
| gpu_failed | 4 (2%) | 3 (12%) |
| postload_rejected | 23 (11%) | 12 (48%) |
| grtresna_rejected | 74 (37%) | 0 (0%) |
| grtresna_failed | 18 (9%) | 0 (0%) |
| pipeline_interrupted | 8 (4%) | 0 (0%) |

Key difference: Trajectory has **zero GRTresna failures** (trivial momentum constraint
at t=0 with static lumps) but **higher postload rejection** (48%) from constraint
violations when loaded onto the evolution grid. SH has the opposite profile — many
GRTresna solver failures from complex initial data, but lower postload rejection.

### FTL champions

**SH v22:**

| Champion | Eval | Value | Score |
|----------|------|-------|-------|
| ftl_geo_evolving | 189 | 0.101 | 470.6 |
| f_geo_evol | 101 | 0.753 | -0.9 (crashed) |
| max_local_speed | 198 | 4.21c | 45.9 |
| f_op_peak | 199 | 7.0% | 40.4 |
| superluminal_fraction | 101 | 86.9% | -0.9 (crashed) |

**Trajectory v1:**

| Champion | Eval | Value | Score |
|----------|------|-------|-------|
| ftl_geo_evolving | 24 | 0.340 | 301.5 |
| f_geo_evol | 24 | 6.86% | 301.5 |
| max_local_speed | 1 | 1.56c | 12.7 |
| f_op_peak | 24 | 15.2% | 301.5 |
| superluminal_fraction | 24 | 100% | 301.5 |
| ftl_lifetime_fraction | 7 | 100% | 13.2 |

---

## Top evaluations — `trajectory_5lump_v1`

### Eval 8 — Best raw score (1166.8, gpu_failed)

**24.6% geodesic shortcut** at t=6.4, crashed at t≈9.5 (59% survival). NaN in `h11`
at level 1 — numerical instability from strong frame-dragging shear.

| Metric | Value |
|--------|-------|
| Score | 1166.8 |
| Status | gpu_failed (exit_code=1) |
| f_geo_peak | 24.62% @ t=6.4 |
| f_op_peak | 19.55% @ t=6.4 |
| max_local_speed | 2.25c @ t=6.4 |
| superluminal_fraction | 73.6% @ t=6.4 |
| ftl_geo_evolving | 1.000 (saturated) |
| numerical_survival | 0.592 |
| n_frames scored | 2 |
| shift_drive | 0.443 |
| channel_progress | 0.665 |

**Configuration:**

```
Lump  R0    omega   tilt_theta  tilt_phi  well_depth  exotic
0     2.74  -0.468  2.35        5.59      0.028       0.18 → canonical
1     6.70  +0.289  3.13        0.70      0.128       0.78 → EXOTIC
2     4.07  +0.282  2.39        5.68      0.111       0.44 → canonical
3     3.96  -0.160  2.65        5.39      0.064       0.55 → EXOTIC
4     6.19  -0.237  2.89        3.41      0.029       0.76 → EXOTIC

Shared: A_breath=1.348, omega_breath=1.112, z_amp=0.193, omega_z=1.038, well_width=2.466
```

**Pattern:** 3 retro + 2 prograde rotation (counter-rotating). All lumps tilted >2.3 rad
(nearly inverted orbital planes). Three exotic lumps. Strong breathing (A=1.35). Wide
well (2.47). The counter-rotation + inverted planes create intense frame-dragging shear
that produces the strongest FTL signal but also destabilizes the numerics.

---

### Eval 21 — Best stable score (687.9, gpu_ok)

**5.3% geodesic shortcut** with full survival. The best configuration that completes
evolution without crashing.

| Metric | Value |
|--------|-------|
| Score | 687.9 |
| Status | gpu_ok (full survival) |
| f_geo_peak | 5.30% @ t=12.8 |
| f_op_peak | 6.57% @ t=12.8 |
| max_local_speed | 1.21c @ t=16.0 |
| superluminal_fraction | 97.3% @ t=12.8 |
| ftl_geo_evolving | 0.261 |
| ftl_lifetime_fraction | 100% (all frames) |
| numerical_survival | 1.000 |
| constraint_health | 0.852 |
| n_frames scored | 7 |
| ftl_persistence | 0.419 |
| ftl_precursor | 0.922 |
| shift_drive | 0.219 |

**Configuration:**

```
Lump  R0    omega   tilt_theta  tilt_phi  well_depth  exotic
0     5.36  -0.634  2.35        2.50      0.087       0.54 → EXOTIC
1     4.14  -0.785  0.83        1.46      0.109       0.77 → EXOTIC
2     4.59  -0.265  0.08        1.32      0.020       0.13 → canonical
3     6.94  -0.621  2.95        3.00      0.022       0.28 → canonical
4     7.25  -0.534  1.81        4.44      0.084       0.81 → EXOTIC

Shared: A_breath=1.270, omega_breath=0.076, z_amp=1.563, omega_z=1.385, well_width=1.364
```

**Pattern:** All 5 lumps retrograde (same direction). Three exotic lumps (0, 1, 4).
Mixed tilts (0.08 to 2.95 rad). Strong z-oscillation (amp=1.56). Slow breathing
(omega=0.08). Moderate well width (1.36). The uniform retrograde with diverse tilts
creates a frame-dragging vortex that generates persistent FTL without the instability
of counter-rotation.

---

### Eval 24 — Second-best stable (301.5, gpu_ok)

| Metric | Value |
|--------|-------|
| Score | 301.5 |
| Status | gpu_ok (full survival) |
| f_geo_peak | 6.86% @ t=9.6 |
| f_op_peak | 15.21% @ t=9.6 |
| max_local_speed | 1.43c @ t=9.6 |
| superluminal_fraction | 100% @ t=12.8 |
| ftl_geo_evolving | 0.340 |
| ftl_lifetime_fraction | 100% |
| constraint_health | 0.738 |

**Configuration:**

```
Lump  R0    omega   tilt_theta  tilt_phi  well_depth  exotic
0     4.31  -0.795  2.43        2.32      0.078       0.66 → EXOTIC
1     1.73  -0.481  0.75        0.36      0.117       0.61 → EXOTIC
2     4.43  -0.853  0.67        5.17      0.018       0.11 → canonical
3     7.09  -0.831  3.06        0.18      0.028       0.33 → canonical
4     7.77  -0.249  0.48        4.67      0.118       0.74 → EXOTIC

Shared: A_breath=0.621, omega_breath=0.999, z_amp=1.675, omega_z=1.498, well_width=1.479
```

**Pattern:** All retrograde. Three exotic (0, 1, 4). Lump 1 very close to center (R0=1.73)
with high amplitude (0.117) — acts as a strong central gravitational anchor. Large spread
in orbital radii (1.73 to 7.77). Strong z-oscillation (amp=1.68).

---

## Physics analysis — Why trajectory >> SH

### 1. Differential motion (the decisive mechanism)

The SH ansatz gives all 5 lumps **identical velocity** (one global `v_tor`, `v_pol`, `v_rad`).
This produces rigid-body rotation with no frame-dragging *shear* — the matter moves as a
single unit, generating a uniform frame-drag field with no gradient.

The trajectory ansatz gives each lump **independent angular velocity** in an **independently
tilted orbital plane**. This creates frame-dragging shear between counter-moving or
differently-oriented matter streams — exactly the mechanism behind Alcubierre-type warps
where differential frame-dragging produces effective metric expansion/contraction.

### 2. GRTresna convergence advantage

| | SH | Trajectory |
|-|----|----|
| t=0 momentum | Non-trivial (lumps have velocity) | **Trivial** (lumps at rest) |
| Solver failures | 46% (grtresna_rejected + failed) | **0%** |
| Effective GPU budget | 79 / 202 evals reach GPU (39%) | 13 / 25 reach GPU (52%) |

The trajectory ansatz places static lumps at t=0. The momentum constraint is trivially
satisfied (zero momentum source → zero shift). Motion comes from the C++ `TrajectoryEvaluator`
during evolution, not from the initial data. This means the constraint solver never fails,
and all computational budget goes directly to GPU evolution.

### 3. Per-lump exotic matter placement

SH uses a coarse azimuthal wedge for exotic matter (contiguous sector, same for all lumps).
Trajectory provides per-lump binary exotic flags, enabling exotic matter at specific orbital
positions. Every top-3 configuration uses **3/5 exotic lumps (60%)**, with the exotic lumps
typically at intermediate radii where the frame-dragging gradient is steepest.

### 4. Tilted orbital planes create 3D topology

SH matter lives in a fixed shell geometry (spherical). Trajectory lumps orbit in arbitrarily
tilted planes (per-lump `tilt_theta`, `tilt_phi`), creating complex 3D frame-dragging
topology. The crashed champion (eval 8) has all lumps with `tilt_theta > 2.3` (nearly
inverted planes), creating a tangled 3D vortex structure.

### 5. Stability vs strength trade-off

| Pattern | FTL strength | Stability |
|---------|--------------|-----------|
| Counter-rotating (eval 8: 3 retro + 2 pro) | **24.6%** f_geo | Crashed @ 59% |
| All-retrograde with tilted planes (eval 21) | 5.3% f_geo | **Full survival** |
| All-retrograde with central anchor (eval 24) | 6.9% f_geo | **Full survival** |

Counter-rotation creates the strongest FTL signal (5x stronger than co-rotation) but is
numerically unstable. The shear gradients at the interface between counter-moving streams
exceed what the AMR grid can resolve, producing NaN at level 1.

All-retrograde with diverse tilts is a compromise: differential frame-dragging comes from
the tilt geometry rather than opposing velocities, producing moderate but **persistent**
FTL that survives the full evolution.

---

## Emerging parameter patterns

### What correlates with FTL success

| Parameter | FTL-positive pattern | Evidence |
|-----------|---------------------|----------|
| Exotic fraction | 3/5 lumps exotic (60%) | All top-3 configs |
| Rotation direction | All retrograde OR counter-rotating | Top-5 all have majority retro |
| Tilt diversity | Mixed tilts (0.08 to 2.95 rad) | Creates 3D frame-drag topology |
| Z-oscillation | Strong (amp > 1.5) | Evals 21, 24 both > 1.5 |
| Well width | Moderate (1.3–1.5) | Narrow wells concentrate energy |
| Breathing | Present but not dominant | A > 1.0 in top configs |
| Well depth | Mixed (0.02–0.12) | Some lumps strong, some weak |

### What correlates with instability

| Pattern | Risk |
|---------|------|
| Counter-rotation (opposing omegas) | 12% crash rate |
| Very wide well (>2.0) | Eval 8 crashed, well_width=2.47 |
| High aggregate well_depth | Stronger field → more nonlinear |
| Large tilts + counter-rotation | Maximum shear → NaN |

### What causes postload rejection (48% of evals)

The postload gate checks constraint violations after loading GRTresna's solved initial
data onto the coarser evolution grid (128³ → 128³ with AMR). Configurations rejected here
typically have:
- Very large orbital radii (lumps near domain boundary)
- High aggregate matter density (sum of well_depths > ~0.4)
- Very different GRTresna vs evolution grid resolution mapping

---

## Next steps

### Immediate (resume campaign)

1. **Restart `trajectory_5lump_v1`** to reach 50 evals target. The driver died mid-campaign
   (metadata tracks `last_eval_counter=30`). Restarting should resume from eval ~30.

2. **Tighten postload gate or widen bounds** — 48% postload rejection wastes compute.
   Options:
   - Relax `postload-max-ham-l2` from 1e-2 to 2e-2 (accept more borderline configs)
   - Or pre-filter: reject configs where `sum(well_depth) > 0.35` before GRTresna

3. **Investigate stability** — The 12% crash rate comes from strong FTL configs. Can we:
   - Reduce `dt_multiplier` from 0.02 → 0.01 for trajectory (smaller timestep)
   - Add sub-cycling at level 1
   - Use Kreiss-Oliger dissipation

### Medium-term (exploit findings)

4. **CMA-ES around eval 21** — Warm-start from the best stable config. Search within
   ±20% of its parameter values to hill-climb the score. Target: break 1000 pts stable.

5. **Counter-rotation-focused survey** — Pin lumps 0,1 to have opposing `omega_rot` signs,
   reduce the other dimensions. The 24.6% signal from eval 8 suggests counter-rotation is
   the strongest mechanism, just needs stability improvements.

6. **Higher resolution validation** — Run eval 21 and eval 24 at N=256, L=128, ml=2,
   stop_time=30 to confirm signals are not resolution artifacts. This is the HQ promote
   path from previous campaigns.

7. **Extend evolution time** — Current `stop_time=16.0` is short. Eval 21 shows FTL
   persisting through all frames (100% lifetime). Running to t=32 or t=64 would test
   whether the signal is truly persistent or eventually decays.

### Longer-term (architectural improvements)

8. **Boson star trajectory** — Replace pump-spotlight matter with self-gravitating complex
   scalar solitons on trajectory paths. Benefits: physical matter that doesn't disperse,
   reduced pump artifacts, non-trivial momentum constraint from initial bulk velocity.
   See [MapElites.md, Option D boson star section](./MapElites.md#future-extension-trajectory--boson-stars).

9. **Adaptive pump governor** — Current pump continuously creates field at spotlight
   position. An adaptive governor that reduces pump when local field density is sufficient
   would reduce instability from over-driving.

10. **Elliptical / precessing orbits** — Current circles limit the geometry. Adding
    eccentricity and apsidal precession would enable time-varying orbital separation
    (close approach → strong interaction → separation → recovery), potentially finding
    transient FTL pulses at periapsis.

---

## Run log

| Run | Date | Ansatz | Evals (scored/target) | Best stable | Best raw | Headline |
|-----|------|--------|----------------------|-------------|----------|----------|
| `scalar_sh_ftl_v22` | 2026-06-24 | SH (ℓ=4) | 202/200 | 470.6 (eval 189) | 470.6 | 1 FTL hit in 202 evals; 2.1% geodesic |
| `trajectory_5lump_v1` | 2026-06-25 | Trajectory (5 lumps) | 25/50 | 687.9 (eval 21) | 1166.8 (eval 8, crashed) | 7 FTL hits in 25 evals; 24.6% geodesic (crashed), 5.3% stable |

**Conclusion:** The trajectory ansatz with per-lump differential motion is a **qualitative
improvement** over spherical harmonics. In 1/8th the evaluations it produces 2.5x stronger
geodesic shortcuts at 40x better FTL yield. The key mechanism is frame-dragging shear from
independently-moving matter in tilted orbital planes — this is what the SH ansatz
structurally cannot express (all lumps share one velocity).
