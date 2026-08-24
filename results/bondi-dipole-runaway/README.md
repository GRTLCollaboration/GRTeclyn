# Bondi dipole runaway — publishable results pack

A positive-active-mass and a negative-active-mass soliton, released at rest in
full 3+1 numerical relativity with dynamical, constraint-solved matter,
**self-accelerate in the same direction** while both same-sector controls stay
put. This pack holds every light artefact behind that claim: per-cell time
series, the dressed-star tables, the solve/evolution parameters, the code
patches, curated frames and movies, and the derived article tables.

Heavy artefacts (plotfiles, gridinit, HDF5 metric stacks, the movies, and the
full ~250-frame-per-field series) stay in the gitignored `runs/` tree on the
machine that produced them. Cells run with frames keep one still per field every
`dt = 10` here. Everything in this pack is small enough to live in git.

## Headline result

A mixed pair — one positive-mass canonical star, one negative-mass phantom star
of matched |ADM| — accelerates as a unit, and **the effect survives grid
refinement**:

| grid | cells | drift at `t = 200` | separation | acceleration |
|---|---|---|---|---|
| coarse | 128³ | `+2.8815` | `10.000 → 10.003` | `1.463e-04` |
| middle | 192³ | `+3.0139` | `10.000 → 9.930` | `1.611e-04` |
| **fine** | **256³** | **`+3.0016`** | `10.000 → 9.915` | **`1.596e-04`** |

The two finest grids agree to `0.4%` on drift and `0.9%` on acceleration. That
is the shape of a converging quantity: refining does not push the effect toward
zero, it settles it onto a value. Quote the `256³` number.

Independent checks bracket it:

| check | cell | result |
|---|---|---|
| does the box matter? | `L = 128`, doubled | drift within `4%` of `L = 64` |
| does it keep accelerating? | `t = 400` | `a` flat at `1.44e-04` through the final quarter |
| does the pull follow the source? | phantom `0.7995 ×` mass | canonical pull ratio `0.809` vs `0.7995` predicted |
| does mesh refinement change it? | `max_level = 1` | level 1 never tagged; drift matches to `0.001%` |

**And the sign flips with the mass ordering.** Four cells spanning mass ratios
`0.597` to `1.333`:

| `\|M−\|` ÷ `M+` | `1.333` | `1.000` | `0.799` | `0.597` |
|---|---|---|---|---|
| separation, 0 → 200 | 10.00 → **10.60** | 10.00 → 10.00 | 10.00 → 9.41 | 10.00 → 8.62 |
| | **opens** | flat | closes | closes more |

The gap opens, sits flat, or closes strictly according to which star is heavier,
crossing zero at the equal-mass point. That is the test that separates a real
gravitational effect from an artefact: no artefact story predicts a sign reversal
that tracks the mass ordering, still less one that crosses exactly where the
masses match.

And the separation law across four cells at `d = 8/10/12/16` gives
`a ∝ d^−2.051` against `−2` exact, with `a·d²` returning the star's mass to
better than `3.5%` at every separation.

**Two honest limits, stated up front.**

*The same-sign pairs merge, and their per-star story is only measured where
frames exist.* Two positive stars attract: the PP N=256 frames show the two
wells closing `8.75 → merged` by `t ≈ 35`, and the ×7 box-activity growth that
follows is merger ejecta draining out through the sponge — not a boundary
artefact (net inward boundary flux stays at round-off all run; an earlier
version of this pack claimed a boundary wave at `t ≈ 32`, and the flux
measurement retires that). The null verdict is *strengthened*: the pair centroid
moves ≤ `7e-04` over the full `t = 200`, through the merger. What is not yet
measured is the MM side — same ×7 growth, but both stars share one tracking
sector, so whether two negative masses repelled as Bondi predicts awaits the
frames re-run. The *mirror* cell, zero net mass, remains the cleanest null.

*No wave zone is reached.* With shells at `R = 16/24/32/40`, `ψ₄(l=2)` falls as
`r^−4.0`; radiation would give `r^−1`. Every shell is in the near zone, so this
pack reports no gravitational radiation rather than a measurement of it — which
is consistent with a pair whose total momentum is zero and whose mass dipole is
therefore static.

Per-cell numbers, the caveats in full, and how to read a cell name are in
[`campaign/README.md`](campaign/README.md).

## What is where

| path | contents |
|---|---|
| [`MATTER_MODEL.md`](MATTER_MODEL.md) | the bicomplex model, where the sign flip lives, dressed-star initial data, code map |
| [`LAUNCH.md`](LAUNCH.md) | exact launch command and full configuration for every cell |
| `campaign/<cell>/` | every run — the resolution ladder, the separation series, the nulls, the wave-zone box and the long run (see the data dictionary below, and `campaign/README.md` for the guide) |
| `stars/` | dressed-star profile tables `r φ₀(r) α(r)` + the M(ω) family scan |
| `analysis/` | derived tables and the scripts that regenerate them |
| `campaign/<cell>/frames/` | one still per field every `dt = 10`, for the cells run with frames. The movies themselves stay in the gitignored run tree |
| `patches/` | the matter-model modifications this campaign required |

### Data dictionary — `campaign/<cell>/`

Every `.dat` stream carries its own `#` header line naming each column.

| file | what it is | key columns |
|---|---|---|
| `sector_barycenters.dat` | **the trajectory record** — per-sector integrals | 1 `t`, 2 `total_canon`, 3 `bary_x_canon`, 6 `rms_canon`, 7 `total_phantom`, 8 `bary_x_phantom`, 11 `rms_phantom` |
| `confinement.dat` | core health / grip | 3 `peak_activity`, 5 `confined_frac`, 18 `min_chi` |
| `psi4_mode_l2m0.dat`, `psi4_mode_l2_all.dat`, `psi4_directional.dat` | radiation extraction; shells are per-cell — `R = 8, 16` on the `L = 64` cells, `R = 16, 24, 32, 40` on the wave-zone cell | — |
| `shell_profiles.dat` | metric on extraction shells (χ, lapse, K) | — |
| `constraint_norms.dat` | **constraint violation during the evolution** (downsampled to Δt = 0.5; column names are supplied by the pack script, the raw stream has none) | 1 `t`, 2 `L2_Ham`, 3 `L2_Mom`, 16 `Linf_Ham_amr` (worst single point). Columns 7–8 are relative measures — see the caution in the constraint section |
| `energy_conditions.dat`, `curvature_invariants.dat` | NEC/WEC monitors and invariants (downsampled; column names supplied by the pack script) | — |
| `collapse_diagnostics.dat` | **horizon watch** — how deep the gravity well gets and whether a trapped surface appears (downsampled; column names are supplied by the pack script, the raw stream has none) | 1 `t`, 2 `min_lapse`, 3 `min_chi`, 4 `max_abs_K`, 8 `max_ah_r`, 9 `min_theta_plus` |
| `boundary_flux.dat`, `areal_radius.dat`, `ftl_timeseries.dat` | outflow, minimal areal radius, geodesic diagnostics | — |
| `grtresna_params.txt` | the constraint-solve input (couplings, per-lump seeds, tolerances) | — |
| `evolution_params.txt` | the GRTeclyn evolution input (grid, AMR, boundaries, dt) | — |
| `Ham_and_Mom_errors.txt` | per-iteration constraint residuals (last row = final) | — |
| `metadata.json` | provenance: git commit, overrides, solve convergence | — |
| `initial_data.matter.json` | the matter configuration actually seeded | — |

Grid geometry is **not** the same for every cell — it is what the cell name
records. Most cells are `L = 64`, `N = 128`, centre at `32 32 32`, so a star
seeded at `x = ±5` reads as `bary_x = 37 / 27`; the ladder rungs raise `N` to
192 and 256 at the same `L`, and the wave-zone cell doubles `L` to 128 with
`N = 256` (centre `64 64 64`). Barycentre coordinates in the streams are
absolute; drifts quoted anywhere in this pack are `x(t) − x(0)`.

### `analysis/`

| file | contents |
|---|---|
| `summary.csv`, `summary.md` | one row per cell: birth checks → final state |
| `trajectories.csv` | drift / separation / core series for every cell, sampled every Δt = 4 |
| `newtonian_reference.csv` | point-mass Bondi integration with the measured ADM masses, aligned to the NR output |
| `separation_scaling.csv`, `separation_scaling.py` | the separation series: measured drift against a point-mass integration of the *same* configuration, at separations 8 / 12 / 16 — the test that isolates the overlap |
| `convergence_check.csv`, `convergence_check.md` | drift / gap / control artefact across grid resolutions, with the spread between grids |
| `wave_check.csv`, `wave_check.md` | gravitational-wave amplitude on each extraction shell in retarded time, and whether the shells agree |
| `constraint_check.csv`, `constraint_check.md` | constraint violation at both stages — the initial-data solve and the evolution — plus whether refining the grid reduces it |
| `make_tables.py` | regenerates `summary.*` and `trajectories.csv` for the Δx = 0.5 production cells |
| `newtonian_reference.py` | regenerates the point-mass reference (pure stdlib RK4) |
| `convergence_check.py` | regenerates `convergence_check.*` from `campaign/` |
| `wave_check.py` | regenerates `wave_check.*` from `campaign/` |
| `constraint_check.py` | regenerates `constraint_check.*` from `campaign/` |
| `star_family_scan.py` | regenerates `stars/star_family.csv` — the M(ω) table the matched pairing was chosen from (needs the wrapper venv) |
| `star_family_massratio_scan.py` | regenerates `stars/star_family_massratio.csv` — both branches out to ω = 1.0, and the phantom branch's floor on how light its stars get (needs the wrapper venv) |

## campaign/ — the cells behind the result (complete)

Twenty-four cells, one folder each, all on **uniform** grids: the convergence
argument needs a single cell size everywhere, so nothing here uses mesh
refinement — with one deliberate exception, `amrcheck_*`, which switches it on
once to show the choice costs nothing. Each folder carries the same streams and provenance files — the
per-sector time series, the four evolution diagnostics downsampled to
`dt = 0.5`, the elliptic residual history, the solve and evolution parameter
files, and the exact launch environment.

`campaign/README.md` is the detailed guide: how to read a cell name, what every
file is, the per-cell numbers, and the caveats. The short version:

| group | cells | what it settles |
|---|---|---|
| `runaway_pair_d{08,10,12,16}_L64_N128` | 4 | the separation law, `a ∝ d^−2.051` |
| `runaway_pair_d10_L64_N{192,256}` | 2 | the resolution ladder's finer rungs — `d10_N128` above is its base. **The headline** |
| `control_lone_{canonical,phantom}` | 2 | a lone star does not drift |
| `control_pair_{pp,mm}_d10_L64_N128` | 2 | same-sign pairs merge yet their centroid never moves — no runaway |
| `control_pair_pp_d10_L64_N{192,256}` | 2 | the null's own resolution ladder |
| `control_pair_mm_d10_L64_N192` | 1 | the mm null's second rung |
| `control_mirror_mp_d10_L64_N128` | 1 | flipping the pair reverses the runaway exactly |
| `massscale_*` + `massratio_*` (2) | 3 | the pull follows the partner's mass — and reverses with it |
| `amrcheck_pair_d10_L64_N128_lev1` | 1 | mesh refinement changes nothing (six-digit match) |
| `wavezone_pair_d10_L128_N256` | 1 | doubled box, four extraction shells |
| `longrun_pair_d10_t400_L64_N128` | 1 | the acceleration is steady to `t = 400` |
| `stability_canonical_w{075,080,085,090}` | 4 | the star family is stable at all |

### What the campaign shows

1. **The runaway converges.** N=192 and N=256 agree to `0.4%` on drift. A
   discretisation artefact shrinks toward zero under refinement; this does not.
2. **It obeys an inverse-square law.** Across `d = 8/10/12/16`, `a·d²` returns
   the star's mass to better than `3.5%`, tightening monotonically with
   separation exactly where finite-size corrections predict.
3. **It reverses when the pair is mirrored,** to within `0.002%` — a sign flip
   no artefact story reproduces.
4. **It scales with the partner's mass.** Lightening the phantom to `0.7995` of
   matched cuts the canonical's pull to `0.809` of its matched value and leaves
   the phantom's own pull within `2.7%` of unchanged.
5. **It is steady, not explosive.** At `t = 400` the acceleration is still
   `1.44e-04`, flat.

### Two limits this pack does not paper over

**The same-sign cells have no per-star tracking.** The sector splitter assigns
matter by field sign, so both stars of a same-sign pair land in one sector: the
tracker reports a single core at the pair midpoint and `separation = nan`.
Their sevenfold activity growth (onset `t ≈ 40`, peak `t ≈ 95`, flat across
N = 128/192/256) is, where frames allow it to be measured, the *merger* of the
two stars — the PP wells close `8.75 → merged` at `t ≈ 35` and the ejecta later
drains through the sponge, with net inward boundary flux at round-off
throughout. (An earlier revision blamed a boundary wave arriving at `t ≈ 32`;
the flux measurement retires that.) The centroid null holds over the full run:
≤ `7e-04` against the runaway's `+2.88`. The MM pair's per-star fate — repulsion
is the Bondi prediction for two negative masses — is unmeasured pending the
frames re-run, and the flux diagnostic reads zero in a phantom-only box, so it
cannot arbitrate. The mirror cell, zero net mass, is the cleanest null.

**No wave zone is reached.** `r·ψ₄` should be flat in the radiation zone. Across
`R = 16/24/32/40` it falls by a factor of `16`, a power law of `r^−4.0` against
`r^−1` for radiation. Every shell available in this box is near-field. This pack
therefore reports *no measurable gravitational radiation* out to `R = 40`,
which is consistent with a pair of zero total momentum whose mass dipole is
static — not a measurement of a wave.

### Frames and movies

Cells run with frames keep one still per field every `dt = 10` under
`<cell>/frames/`, with a `FRAMES.md` listing the times kept. The full series
(~250 frames per field) and the stitched movies stay in the gitignored run tree
— 19 movies per cell, redrawn against a single fixed colour scale measured over
the whole run so a colour means the same value in every frame.

The `t = 400` cell keeps no stills: it exists for one moving picture, is twice
as long as every other cell, and no number in the analysis is read off it.

### An earlier campaign lived here

Twenty `convA_*`/`boxC_*` cells occupied this folder until 2026-08-22. They
were run before the initial-data grid alignment was fixed, so every star was
born displaced from the centre of its own gravitational well by a fraction of a
cell — an artefact the same size as the signal. They are in git history and
nothing here depends on them.

## Reproducing

Regenerate this whole pack from the run tree:

```bash
bash research/bondi_dipole/pack_runaway.sh          # every cell -> campaign/ + stars/
FRAME_DT=20 bash research/bondi_dipole/pack_runaway.sh   # fewer stills
```

Safe to re-run at any time, including while cells are still evolving: each cell
folder is rebuilt from scratch, and cells with no time series yet are skipped
with a note. It copies only small artefacts and scrubs absolute machine paths at
runtime (`grteclyn-wrapper/src/grteclyn_wrapper/packaging/scrub_paths.py`) — no
host, user, or site identity enters git.

`FRAME_DT` sets the frame spacing (default 10) and `FRAMES_SKIP` names cells
whose stills are not packed at all; the `t = 400` cell is skipped by default.
Frame thinning itself lives in `research/bondi_dipole/thin_frames.py`.

Two earlier packers, `pack_results.sh` and `pack_campaign.sh`, are **superseded
and refuse to run** — they read run trees that no longer exist. They are kept
for the reasoning in their headers.

To re-run the physics itself, see [`LAUNCH.md`](LAUNCH.md).

## Provenance

- Initial data: **GRTresna** (CTTKHybrid, `BosonStarBH` example, complex scalar
  matter with per-lump signs), 32 MPI ranks on CPU, one solve at a time.
- Evolution: **GRTeclyn** (`RadialRecipe`, CCZ4 + bicomplex scalar matter),
  single rank, one GPU per cell, no mesh refinement.
- Working notes (narrative, not part of this pack): [`research/bondi_dipole/docs/`](../../research/bondi_dipole/docs/)
  (working copy: `research/bondi_dipole_debug.md`).
- Code state: wrapper commit recorded per cell in `campaign/<cell>/metadata.json`;
  the GRTresna matter modifications are in [`patches/`](patches/).
