# Pre-GPU pipeline — findings log

Status 2026-07-29. Branch `feature/interstellar`.
Scope: everything **upstream** of GPU evolution — MAP-Elites/CMA-ES proposal,
genome→config decode, GRTresna elliptic solve, Chombo HDF5→`.gridinit`
conversion, matter metadata, and the acceptance gates.

Companion to [`Debug.md`](Debug.md), which covers the **post-GPU** stage
(diagnostics, scoring, probes). That log found 22 bugs. This one exists because
the same audit had never been run on the stage that *produces the candidates*.

**Fix pass complete.** Nine of the ten findings are fixed with a regression test
that fails before the change and passes after; PG-3 is partly fixed and its
remaining work is a measurement, not a patch. Three things were deliberately not
done and are named where they arise: the composite constraint threshold (no data
to derive it from yet), whether the gate run should refine (costs GPU time per
candidate), and whether the 5-lump cap is still right (widens the state vector).

**Standing rules.**
1. Nothing here has been written into `research.tex`.
2. CONFIRMED means both sides of the interface were read in code **and** the
   effect was reproduced in stored run data. PLAUSIBLE means mechanism traced,
   effect not yet measured.
3. A claim about what the search *found* is only as good as the encode/decode
   round trip. A claim about what the geometry *is* survives a broken round trip
   — keep the two separate when triaging.

## Severity scale

| level | meaning |
|---|---|
| **S1 CRITICAL** | invalidates a published claim, or corrupts the search process itself |
| **S2 MAJOR** | silently wrong physics, or a search dimension that does nothing; measured numbers may survive, the narrative does not |
| **S3 MINOR** | real defect, bounded and quantifiable impact |
| **S4 hygiene** | no effect on any result; fix when convenient |

## Findings at a glance

| id | severity | status | one line | test |
|---|---|---|---|---|
| PG-1 | S1 | **FIXED** | trajectory genome is re-derived from decoded physics each generation; rotation/drift not heritable | `test_genome_roundtrip.py` |
| PG-2 | S2 | **FIXED** | "dynamical" boson-shell mode searches no velocity, hides a fixed 0.35c current | `test_shell_static_coherence.py` |
| PG-3 | S1 | **PARTLY FIXED** — gate reads the composite columns; **threshold still uncalibrated, and no gate run has ever refined** | no instrument has measured constraint satisfaction where the matter is | `test_postload_gate_composite_norms.py` |
| PG-4 | S3 | **RESOLVED** | paper's 3e-2 is correct; the 1e-2/3e-2/0.1 spread across campaigns is the real issue | — |
| PG-5 | S3 | **FIXED** | RNG stream duplication (22 workers, 18 streams) and resume replay | `test_pipeline_hygiene.py` |
| PG-6 | S3 | **FIXED** | exotic-injection knob is a silent no-op for trajectory campaigns | `test_pipeline_hygiene.py` |
| PG-7 | S4 | **FIXED** | retrograde fold makes the omega axis non-identifiable | `test_retrograde_axis.py` |
| PG-8 | S3 | **FIXED** | NaN descriptor crashes the completion callback, leaks a pipeline slot | `test_pipeline_hygiene.py` |
| PG-9 | S1 | **FIXED** (parts 1–2; part 3 is a physics decision) | 8-lump campaign declared 56 trajectory dimensions; GRTeclyn silently drives 5 lumps, 21 dimensions inert | `test_param_contract.py -k eight_lump` |
| PG-10 | S2 | **FIXED** | `.gridinit` was loaded by column position on the C++ side and by column name on the Python side | `test_cross_code_contract.py` |

Executable checks for all of the above are described in **PART F**.

**Resume replay (part of PG-5) is NOT fixed.** Every campaign still starts from
the same fixed seed 7, so a resumed run replays the original run's early noise.
That is a campaign-script decision, not a code defect.

---

# PART A — CONFIRMED FINDINGS

## PG-1 [S1 CRITICAL, CONFIRMED] — trajectory genome round-trip asymmetry: rotation and drift are not heritable

**The single most consequential pre-GPU defect found.** It disables 16 of the
39 search dimensions and it wrote a false sentence into the paper.

### Mechanism

The genome coordinate `trajectory_lump<k>_omega_rot` ∈ [−1, 1] is a
**normalized fraction of the speed budget**. Decode converts it to a physical
angular velocity:

```
omega_phys = omega_norm · v_max / R0        # v_max = 0.3 default, R0 ≈ 2.6–5.6
v_rad_phys = v_rad_norm · v_max
```

`_clamp_trajectory_speed` writes the **physical** result back into the same
dict key the genome coordinate came from
([`optimize/candidates.py:171`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/optimize/candidates.py#L171)),
and that dict is what gets stored as the elite's genome
([`qd_search/driver.py:353`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/qd_search/driver.py#L353),
`params = {d.param_key: float(overrides[d.param_key]) for d in dims}`).

Every consumer then reads the physical value **as if it were the genome
coordinate**, with no inverse transform anywhere:

| re-encode site | file:line | consequence |
|---|---|---|
| elite mutation (~60% of samples) | [`qd_search/sampling.py:96`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/qd_search/sampling.py#L96) | child rotation = noise, parent contributes ~2% |
| feasible-box bounds | [`qd_search/sampling.py:43-52`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/qd_search/sampling.py#L43-L52) | box collapses onto the contracted range |
| CMA-ES warm start | [`optimize/candidates.py:286-298`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/optimize/candidates.py#L286-L298) → [`:242`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/optimize/candidates.py#L242) | CMA-ES **starts** 13× contracted |
| near-miss / shadow pool | `qd_search/pre_gpu_archive.py:36`, `pre_gpu/near_miss_pool.py:37-44` | same |
| `--seed-eval-dirs` | `qd_search/driver.py:800-803` | same |
| surrogate training | `optimize/eval.py:46-67` | trains on physical X, screens genome X |

The contraction factor is `R0/v_max` ≈ **8.7–18.6× per round trip**, while
mutation noise stays at full scale (σ = 0.15 × range = 0.3). The parent's
rotation signal is therefore 5–27× **below** the noise it is buried in.

The developers half-knew: HQ replay explicitly skips re-conversion
([`scripts/campaigns/hq/replay_eval.py:403-408`](../../grteclyn-wrapper/scripts/campaigns/hq/replay_eval.py#L403-L408)).
The search loop itself was never given the same treatment.

### Evidence — heritability is exactly zero

`runs/grtresna_qd/qball_traj_bicomplex_8lump_v1/trajectory.jsonl`, 200 evals.
Reconstructing the implied genome coordinate `|omega_phys|·R0/v_max` from what
was actually stored:

| evals | implied \|omega_norm\| | implied \|v_rad_norm\| | expected if… |
|---|---|---|---|
| 0–24 (random init) | **0.490** | 0.417 | 0.500 = uniform random ✓ |
| 25–74 | **0.235** | 0.266 | 0.239 = pure mutation noise |
| 75–149 | **0.260** | 0.260 | 0.239 |
| 150–199 | **0.241** | 0.252 | 0.239 |

The signal dies at the exact transition from random seeding to mutation and
never recovers over the remaining 175 evaluations. Uniform init gives
E|x| = 0.500; pure reflected noise at σ = 0.3 gives E|x| = 0.239. The run sits
on the second value from eval 25 onward. **MAP-Elites was doing random search
in these dimensions, not evolution.**

Search-space bounds confirm the axes are normalized (`metadata.json`:
`{'key': 'trajectory_lump0_omega_rot', 'lower': -1.0, 'upper': 1.0}`, all 8
lumps, both axes).

### Evidence — it reached the paper

`research.tex:811` presents candidate 146's near-zero rotation as a discovery:

> "CMA-ES bought that by collapsing onto a quasi-static scaffold---all five
> angular velocities below $5\times10^{-3}$---which keeps the lumps in place
> long enough for the gated terms to pay out."

and `:814`: *"The run is converged rather than merely lucky."*

From the surviving replay metadata
(`runs/grtresna_promote/bcma_rm_L128_N256_t30_hq_eval000146/metadata.json`,
values identical in its `params.txt`, so this is what GRTresna and GRTeclyn
actually received):

| lump | R0 | omega_phys | **implied genome** | v_t = \|ω\|R0 | % of v_max budget | R0/v_max |
|---|---|---|---|---|---|---|
| 0 | 3.69 | −0.001736 | 0.0213 | 0.0064 | 2.13% | 12.3× |
| 1 | 5.50 | −0.000723 | 0.0133 | 0.0040 | 1.33% | 18.3× |
| 2 | 2.83 | −0.003101 | 0.0293 | 0.0088 | 2.93% | 9.4× |
| 3 | 5.59 | −0.000919 | 0.0171 | 0.0051 | 1.71% | 18.6× |
| 4 | 2.61 | −0.004927 | 0.0428 | 0.0128 | 4.28% | 8.7× |

**Mean implied genome coordinate = 0.0248 on a [−1, 1] axis — 2.5% of the
half-range from dead zero.** CMA-ES (σ₀ = 0.05) never explored rotation; it was
*started* at a contracted coordinate and could not escape it. The earlier
retained eval (`sources/qball_traj_bicomplex_cmaes_v1/eval_000063`) sits in the
same band (implied 0.004–0.099), confirming the whole CMA-ES run was confined
to ≲5% of the rotation axis.

**The clincher.** `research.tex:1101-1102` says candidate 146 is *"roughly an
order of magnitude quieter than MAP-Elites champion 87."* The predicted shrink
from exactly one warm-start re-encode is **mean 13.5×**. The paper measured the
bug's contraction factor and wrote it down as a physical finding.

### Consequence

* **"Quasi-static scaffold" must be withdrawn as a result.** It is not what
  CMA-ES found; it is what the re-encode did to the starting point.
* **"The run is converged rather than merely lucky" (`:814`) does not follow.**
  The eight best genomes agreeing to 4 score points is equally explained by ten
  of their dimensions being frozen near zero.
* **The campaign's stated primary mechanism was never searched.** The
  qball_trajectory campaign exists to explore counter-rotating lumps →
  frame-dragging shear (`scripts/campaigns/qball_trajectory/run.sh:22-27`).
  Rotation was re-randomized every generation.
* **What SURVIVES:** the 20.21% headline. HQ replay does not re-convert, so the
  evolved geometry is the one `params.txt` specifies, and the null-geodesic
  measurement of it stands. **The number is real; the account of how it was
  discovered is not.** Same split as the metric-cache bug in `Debug.md` §15 —
  measurement sound, provenance broken.

### Verification (reproduce in ~2 s)

```bash
cd runs/grtresna_qd/qball_traj_bicomplex_8lump_v1
python3 -c "
import json,statistics as st
recs=[json.loads(l) for l in open('trajectory.jsonl') if l.strip()]
rows=[]
for r in recs:
    ov=r.get('overrides') or {}; vmax=float(ov.get('trajectory_v_max',0.3) or 0.3)
    om=[abs(float(ov[f'trajectory_lump{k}_omega_rot']))*float(ov[f'trajectory_lump{k}_R0'])/vmax
        for k in range(8) if ov.get(f'trajectory_lump{k}_omega_rot') is not None and ov.get(f'trajectory_lump{k}_R0')]
    if om: rows.append((r.get('eval'), st.mean(om)))
rows.sort()
for lo,hi,lbl in [(0,25,'init'),(25,75,'25-74'),(75,150,'75-149'),(150,10**9,'150+')]:
    s=[v for e,v in rows if lo<=e<hi]
    if s: print(f'{lbl:8s} n={len(s):3d} implied |omega_norm|={st.mean(s):.3f}')
print('uniform init -> 0.500 ; pure noise -> 0.239')"
```

---

## PG-2 [S2 MAJOR, CONFIRMED (conditional reachability)] — "dynamical" boson-shell mode searches no velocity, and hides a fixed 0.35c current

**Files:** [`cli/grtresna_context.py:58-66`](../../grteclyn-wrapper/src/grteclyn_wrapper/cli/grtresna_context.py#L58-L66),
[`optimize/spaces.py:319-343`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/optimize/spaces.py#L319-L343),
[`optimize/config.py:44-49`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/optimize/config.py#L44-L49).

`shell_static` defaults to **True** unless explicitly pinned to 0. When True,
`spaces.py` **strips the four velocity dimensions from the search space** — but
`grtresna_shell_static` itself remains a searched dimension in [0, 1], and
decode defaults it to **0 = moving**. Any candidate decoding to non-static then
picks up a hard-coded `grtresna_shell_toroidal_velocity = 0.35`, which is
neither searched nor recorded as a search dimension.

`scripts/campaigns/boson_star/ftl_shell_run.sh:37-42` documents `PIN_DIMS=""`
as *"dynamical (velocity/omega searched)"* — with no `shell_static` pin, this
is exactly the broken configuration.

**Failure scenario.** A run advertised as a dynamical boson-shell search
explores zero velocity dimensions, while roughly half its candidates carry an
invisible 0.35c toroidal current absent from every trajectory record. Any
matter-class comparison drawn from that campaign is confounded.

**Not reached by:** default boson campaigns (which pin `shell_static=1`) or any
scalar-shell campaign. Same trap applies to `OBJECTIVE_MODE=critical_collapse`
launched outside `splash/run.sh`.

**Verify:** dry-run `ftl_shell_run.sh` with `PIN_DIMS=""`; `metadata.json`
`search_space` will contain no `*_velocity` key while the emitted `params.txt`
carries velocity 0.35.

**FIXED.** A space that searches no velocity now decodes as static outright,
regardless of the `shell_static` bit, and the flag is dropped from such a space
so no axis is spent on a decision that no longer exists. Confirmed against the
old decoder: the same overrides used to yield `velocity = (0, 0, 0.35)` per lump
and now yield `(0, 0, 0)`.

One thing the obvious version of this fix gets wrong, recorded because it cost a
test: the v1 boson-shell builder puts the currents in its `skip` set and *then*
re-adds boson-specific ones when `static=False`. Deciding on `skip` therefore
strips the flag from moving spaces too. The check has to run on the finished
dimension list.

---

## PG-3 [S1 CRITICAL, mechanism CONFIRMED / magnitude PENDING] — no instrument has ever measured constraint satisfaction where the matter is

Carried over from [`Debug.md` §4 item 0](Debug.md) and confirmed against the
paper's claim inventory. Recorded here because **it is a pre-GPU acceptance
problem**, not a diagnostics problem.

Every "the initial data satisfies the constraints" statement in `research.tex`
traces to exactly two instruments:

1. **GRTresna's own global convergence percentage** — 5%/5% threshold
   ([`search/grtresna_convergence_gate.py`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/grtresna_convergence_gate.py)),
   measured **pre-conversion**, i.e. never on the `.gridinit` that is actually
   evolved.
2. **The post-load gate** ([`projection/postload_gate.py`](../../grteclyn-wrapper/src/grteclyn_wrapper/projection/postload_gate.py)),
   which reads `constraint_norms.dat` **cols 2–3** — the level-0-only,
   unmasked, whole-domain-diluted norm that `Debug.md` §4 item 0 has already
   ruled unquotable as an accuracy figure.

Since ~99.9% of the domain is near-vacuum, the peak `|Ham|` runs ~218× the
composite L2 and the refined-region norm ~30× the domain mean. **A constraint
violation localized exactly at the lumps — the physically dangerous case — can
pass a 1e-2 gate while being 30–200× larger than the number printed.** The
first composite readings (`L2_Ham_amr_rel` ≈ 0.22–0.40 even at t ≤ 2, i.e. the
residual is comparable to the terms it must cancel) indicate this is live, not
hypothetical.

**Paper text at risk:** `:208-211` ("the Hamiltonian and momentum constraints
hold before any evolution" — the strongest ID claim in the paper),
`:181-183` ("constraint-clean data"), `:468-480` (gate values 5.9e-3 / 1.1e-4
and "establishes numerical control and source-model consistency"),
`:1524-1528`, `:1467-1472` (Richardson p_H ≈ 3.3 — certifies convergence of the
level-0 restriction, not of the solution).

**Open:** cols 12–17 (`L2_Ham_amr*`, commit `77af4b02`) exist in the binary but
the gate does not read them, and **no paper-tier run has produced them** —
campaigns A–F all predate the columns and have 11.

### What changed, and the two things that did not

The gate now reads the composite columns and reports all three; a legacy 11-column
file still parses and says so. What the check ordered first — *confirm the gate
run actually has `max_level > 0` and has regridded* — came back **negative**, and
that reshapes the rest of the item.

Every stored `postload_gate` episode (ten of them, across
`qball_traj_bicomplex_8lump_v1` and `qball_traj_bicomplex_cmaes_v1`) logs:

```
grid_places() time: 0.009824038 new finest: 0
STEP = 0 TIME = 0 : REGRID  with lbase = 0
  Level 0   64 grids  2097152 cells  100 % of domain
```

No level 1, ever. Three things combine: `ChiTagger` at `regrid_threshold = 0.01`
does not fire on freshly projected data; `regrid_interval = 16` means the initial
hierarchy build is the only opportunity; and the run is **one coarse step long**
(`L=64`, `N=128` → `dx = 0.5`; `dt = 0.02·dx = 0.01`; `stop_time = 0.01`).

Two consequences:

1. **On today's gate run the composite columns are not composite.** They describe
   the same single level 0 as cols 2–3, differing only by the four dropped
   boundary cells and by reporting the undiluted peak `Linf_Ham_amr`. That peak
   is still worth having — it is the quantity a localized violation moves and the
   domain mean does not — but the ~30× refined-region figure this item was built
   around is not being measured.
2. **Col 17 reads as perfect on exactly the runs where it means nothing.** With
   no level ≥ 1 the refined-region sum is written as `0.0`. Gating on it would
   have passed every candidate unconditionally. `has_refined_region` is False in
   that case and the gate refuses to score it.

**Still open, and deliberately so:** the composite thresholds are unset. There is
no run from which to derive one — the very files that would calibrate it are the
11-column ones. The gate reports the values into `notes` on every run so the data
accumulates; setting the number is a decision to make on that data, not before
it. And making the composite columns measure a real hierarchy needs the gate run
to refine, which costs GPU time per candidate and changes what the gate is — an
operator call, not a code fix.

---

## PG-4 [S3 MINOR, RESOLVED — paper is correct; the spread is the real issue]

`research.tex:474-476` states the bicomplex campaigns used **3×10⁻²** for the
Hamiltonian residual and 10⁻² for momentum.
[`postload_gate.py:22-23`](../../grteclyn-wrapper/src/grteclyn_wrapper/projection/postload_gate.py#L22-L23)
defaults **both** to 1e-2.

**Resolved: the campaigns override the default, and the paper is right.** Both
`qball_trajectory/run.sh:156` and `cmaes_run.sh:93` export
`POSTLOAD_MAX_HAM_L2="${POSTLOAD_MAX_HAM_L2:-3e-2}"`, and `run.sh` never sources
`search_common.sh`, so 3e-2 wins for the bicomplex campaigns. No correction to
the paper is needed on this point.

What the check did turn up is worth keeping: the threshold is **not** uniform
across campaigns. `run_shear.sh:88` and `run_lentz.sh:57` use **0.1** — a 10×
spread, with the code default of 1e-2 as a third value for anything that
overrides nothing. All three are applied to the *same* diluted norm (PG-3), so
the spread is not a considered per-campaign tolerance; it is three unrelated
choices that were never compared. Fold into the PG-3 fix rather than
adjudicating separately: once the norm measures where the matter is, the
thresholds have to be re-derived anyway.

---

## PG-5 [S3 MINOR, CONFIRMED] — RNG stream duplication and resume replay

[`qd_search/driver.py:630-648`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/qd_search/driver.py#L630-L648).
Seeds spawned = `total_slots + max_concurrent + 4` = **18**; pipeline worker
threads = `slots + cpu + scoring_workers` = **22** (`driver.py:711-713`).
Workers 19–22 take `idx % 18` and duplicate child streams 0–3 — two concurrent
workers draw identical mutation noise (candidates correlated, not identical,
since parents differ).

Additionally, resume respawns from the same `SeedSequence(seed)`, restarting
every stream at its beginning. The seed is fixed at **7** by default in both
`qd/run.sh:25` and `cmaes/run.sh:18` for every campaign, so resumed runs replay
the original run's early noise sequence.

**Duplication FIXED; replay NOT.** Each worker now mints its own child stream on
demand, so the count cannot disagree with the worker count again. The replay is a
campaign-script property — the seed is literally written as `7` in both launchers
— and changing it silently would make old runs unreproducible, so it is left for
whoever owns those scripts.

---

## PG-6 [S3 MINOR, CONFIRMED] — exotic-injection knob is a silent no-op for trajectory campaigns

[`optimize/candidates.py:344-371`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/optimize/candidates.py#L344-L371).
`_force_exotic_template` matches only `grtresna_ring_*` / `grtresna_shell_exotic_fraction`
/ `grtresna_lump{k}_*` (via `_LUMP_KEY_RE = grtresna_lump(\d+)_`). Trajectory
campaigns use `trajectory_lump{k}_*`, which never matches.

`cmaes/run.sh:21` sets `EXOTIC_INJECTION_FRACTION=0.1`, so in every qball
CMA-ES run 10% of slots silently degrade to plain jittered warm vectors while
the logs imply exotic templates were injected.

**FIXED.** Adding the pattern to the regex alone would not have been enough: the
template then writes `grtresna_lump{k}_{suffix}` keys a trajectory space does not
contain, so it would still have done nothing while looking fixed. The trajectory
branch sets `trajectory_lump{k}_exotic` and nothing else — the orbit is a
searched dimension, and a template that overwrote it would be deleting search
results rather than seeding them.

---

## PG-7 [S4 hygiene, CONFIRMED] — retrograde fold makes the omega axis non-identifiable

[`optimize/candidates.py:177-206`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/optimize/candidates.py#L177-L206).
`_enforce_retrograde` negates positive omega at decode, so genome `+x` and `−x`
decode identically. CMA-ES covariance and QD mutation therefore waste half the
omega axis, and mutations crossing zero are folded back. Correct given the
stated intent (`trajectory_retrograde_only = 1` is set in the candidate-146
`params.txt`), but it compounds PG-1: stored values are always ≤ 0.

Cleaner: search `|omega|` with the sign pinned.

**FIXED**, by narrowing each omega axis to its non-positive half when
`trajectory_retrograde_only` is set — same identifiability, without adding a
dimension or changing the stored genome layout. The fold stays in the decoder for
pinned and legacy values, but on a retrograde campaign it is now provably a
no-op: decoding a genome from the narrowed box gives the same result with the
fold enabled and disabled, at every point along the axis.

---

## PG-8 [S3 MINOR, PLAUSIBLE] — NaN descriptor crashes the completion callback and leaks a pipeline slot

[`qd_search/descriptors.py:430-431`](../../grteclyn-wrapper/src/grteclyn_wrapper/search/qd_search/descriptors.py#L430-L431):
`_bin_index` calls `math.floor(value*bins)` → `ValueError` on NaN. The gw_beam
axis can pass NaN through (`descriptors.py:352-354`: `max(nan, floor)` → NaN →
`log10` → `np.clip` preserves NaN). The exception propagates inside
`_on_eval_complete` **before** `pipeline.submit(new_x)` (`driver.py:650-706`),
so one in-flight slot is permanently lost per occurrence
(`eval_pipeline.py:183-190` decrements counters but nothing resubmits).

ftl_lifetime / speed / channel axes are NaN-guarded; only Ψ₄-NaN reaches this,
i.e. gw_beam campaigns only. The audited qball run had 0 occurrences in 200
records.

**FIXED**, on both halves. `_bin_index` maps NaN to bin 0 and ±inf to the ends
rather than raising — the raw value stays in `descriptor_details`, so nothing is
disguised as a real measurement. Independently, the resubmit that keeps the
pipeline full moved into a `finally`, so no exception raised while recording a
result can cost a slot. The second half matters more than the first: it closes
the whole class, not the one instance.

---

## PG-9 [S1 CRITICAL, CONFIRMED] — the 8-lump campaign searched 8 lumps; GRTeclyn evolved 5

The `qball_traj_bicomplex_8lump_v1` campaign declares **56 trajectory search
dimensions**: lumps 0–7 × 7 parameters (`R0`, `omega_rot`, `v_rad`, `phase0`,
`tilt_theta`, `tilt_phi`, `exotic`). Its `metadata.json` lists all 56, its
`archive.json` records elites over all 56, and every emitted `params.txt`
carries all 64 keys (8 params × 8 lumps).

GRTeclyn reads five.

[`SimulationParameters.hpp:122-128`](../../Examples/RadialRecipe/SimulationParameters.hpp#L122-L128):

```cpp
int n_traj = 5;
pp.load("trajectory_num_lumps", n_traj, 5);
if (n_traj > GRTRESNA_MAX_INDEPENDENT_SCALARS)
    n_traj = GRTRESNA_MAX_INDEPENDENT_SCALARS;   // = 5, silently
trajectory_params.num_lumps = n_traj;
```

and the per-lump read loop at `:144` runs `k < trajectory_params.num_lumps`,
i.e. `k < 5`. The config says `trajectory_num_lumps = 8`; the clamp takes it to
5 without a warning, a log line, or a nonzero exit. **`trajectory_lump5_*`,
`trajectory_lump6_*` and `trajectory_lump7_*` — 24 keys, 21 of them wired to
live search dimensions — are written into every `params.txt` and read by
nothing.**

**This is not the same as the lumps being absent.** GRTresna is *not* capped: its
`num_lumps = 8` is honoured and `lump0..lump7_*` are all read, so all eight
lumps' matter is genuinely present in the initial data and contributes to the
constraint solve. What is missing downstream is the *trajectory driving* — the
moving pump well — for lumps 5, 6 and 7. They are seeded and then never driven.

### Consequence

Scoped the same way as PG-1: this breaks **provenance**, not measurement.

- **Broken:** any statement about which lump trajectory the search discovered.
  37.5% of the declared search space (21 of 56 dimensions) had zero effect on
  the evolution, so the archive's coverage and the elites' descriptor positions
  are those of a 35-dimensional search reported as a 56-dimensional one.
  QD "coverage" over inert axes is coverage of nothing.
- **Intact:** the measured behaviour of any specific accepted config. The
  evolution did what its `params.txt` specified; it simply specified less than
  the search believed it had specified.

### Fix status

**FIXED, parts 1 and 2.** Both `SimulationParameters.hpp` clamps now
`amrex::Abort` with the requested and available counts; `apply_rl_actions`
asserts. On the wrapper side, `scalar_slot_cap()` reads
`GRTRESNA_MAX_INDEPENDENT_SCALARS` out of the C++ header (currently 5) and both
trajectory search-space builders refuse an over-cap lump count before any GPU
time is spent — including `grtresna_trajectory_boson_search_space`, which is the
builder the 8-lump campaign actually used and which a first pass missed.

**Part 3 is not a code change.** Whether 5 is still the right cap is a physics
and cost decision: raising it widens the state vector by 2 components per extra
lump. Left to the operator.

### Verification

```bash
cd grteclyn-wrapper
PYTHONPATH=src python -m pytest tests/core/test_param_contract.py \
    -k eight_lump -q
```

The checker reports 24 `indexed-key-beyond-cap` findings for that campaign and
none for the 5-lump runs.

### Note on how this was missed twice

The first pass of the params-contract checker **suppressed exactly this
finding**: `trajectory_lump` was on a `_DYNAMIC_KEY_PREFIXES` allowlist, because
those keys are assembled at runtime with an `ostringstream` and the static scan
cannot resolve them. "Assume a reader exists" is the same mistake as a silent
default — it converts *unknown* into *reassuring*. The checker now reads
`GRTRESNA_MAX_INDEPENDENT_SCALARS` out of the C++ header and flags any indexed
key past the cap.

---

# PART B — FIX PLAN, PRIORITY ORDER

### 1. PG-1 — store the genome, stop re-deriving it *(S1, blocks any claim about what the search found)* — **DONE**

Both halves were implemented, not just the cheap one: the raw optimizer vector
is now stored in every trajectory record and in `Elite`, **and** the inverse
transform exists so legacy records without a stored genome still decode
correctly. All seven re-encode sites read the genome first.

One correction to the plan below: the regression test it specifies —
`E(D(x)) == x` — **cannot hold and should not**. The decoder contains three
genuine physical projections (speed disk clamp, spiral-radius floor, retrograde
fold), so for a genome that a cap binds, `E(D(x))` lands on the projected point,
not the original. The tests assert the correct invariant instead:
`E(D(x)) == P(x)` where `P = E∘D`, `P(P(x)) == P(x)`, `P(x) == x` whenever no cap
binds, and the exact clamped value for the one genome that does bind — which is
stricter than the original wording, not weaker.

**Root cause is the `Debug.md` §17 recurring shape:** a value derived across a
transform boundary is read back as truth while **the authoritative source sits
unused** — the raw optimizer vector `x` is right there in `job.payload`.

**Preferred fix (robust):** persist the raw genome vector alongside the decoded
overrides, and make every re-encode site read the genome, never the overrides.

* `qd_search/driver.py:353` — store `record["genome"] = list(x)` and
  `Elite.genome = list(x)` in addition to `params`.
* `sampling.py:_mutate_params`, `_feasible_bounds`,
  `_feasible_bounds_from_params` — take the genome.
* `candidates.py:_load_warm_start_vectors` — read `rec["genome"]` when present;
  fall back to `_overrides_to_vector` **with the inverse transform applied**
  for legacy trajectories.
* `pre_gpu_archive.py`, `near_miss_pool.py`, `optimize/eval.py:_collect_training`
  — same.

**Minimum viable fix:** add the inverse (`omega_norm = omega_phys·R0/v_max`,
`v_rad_norm = v_rad_phys/v_max`) inside `_overrides_to_vector` and
`_mutate_params`. Cheaper, but leaves the two representations sharing one key
name — the exact ambiguity that caused this.

**Regression test that fails before and passes after** (the `Debug.md` §17
seventh-habit rule — *"physics unchanged" is not "defect removed"*):
round-trip a vector `x` with non-zero rotation through
`_vector_to_overrides` → `_overrides_to_vector` and assert
`x_out ≈ x_in` to 1e-12 for every trajectory dimension. **Assert on the
recovered genome, not on physics-neutrality.**

**Then:** re-run one qball QD campaign and confirm the heritability statistic
(PG-1 verification block) stays near 0.5 rather than collapsing to 0.24.

### 2. PG-3 — make the gate measure where the matter is *(S1, gates the paper)* — **PARTLY DONE**

* ~~Point `postload_gate.py` at `constraint_norms.dat` **cols 12–17**~~ **DONE.**
  `read_constraint_metrics` now parses `L2_Ham_amr` (col 12), `Linf_Ham_amr`
  (col 16) and `L2_Ham_amr_ref` (col 17); the gate reports all three and
  enforces whichever of them the caller has configured a threshold for.
* **The threshold is deliberately left unset.** No run has yet produced these
  columns, so any number written today would be invented, not measured. Until
  `PostLoadGateConfig.max_hamiltonian_l2_amr` / `max_hamiltonian_linf_amr` are
  set, the gate's pass/fail is unchanged and the composite values are written
  into `notes` so the next batch of gate runs supplies the calibration data.
  **This is the remaining work on PG-3.**
* **The precondition check failed.** All ten stored `postload_gate` episodes log
  `grid_places() ... new finest: 0` — a single level 0 covering 100% of the
  domain, no level 1 ever built. The gate template asks for `max_level = 2`
  (campaign overrides cut it to 1), but `ChiTagger` at `regrid_threshold = 0.01`
  does not fire on freshly projected data, and `regrid_interval = 16` against a
  run that is **one coarse step long** (`dx = 0.5`, `dt = 0.02·dx = 0.01`,
  `stop_time = 0.01`) means the only chance to refine is the initial hierarchy
  build. So today the composite columns describe the same single-level grid as
  cols 2–3, gaining only the dropped boundary cells and the undiluted peak —
  which is still strictly more than the domain mean, but not what PG-3 wanted.
* **Col 17 is a trap on such a run and is now handled as one.** With no level ≥ 1
  the refined-region norm is written as exactly `0.0`, which reads as flawless
  and means *unmeasured*. `ConstraintMetrics.has_refined_region` is False in that
  case, the gate says so in `notes`, and the value is never treated as a pass.
* **Open decision for the operator:** making the composite columns mean what this
  item wants requires the gate run to actually refine — lower `regrid_threshold`
  for the gate, or run it long enough to regrid. Both cost GPU time per candidate
  and change what the gate measures, so neither was done unilaterally.
* Leave the governor's input (cols 2–3) untouched — changing it is a physics
  change, per `Debug.md` §4 item 0. **Untouched.**

### 3. PG-4 — resolve the threshold discrepancy *(S3, but blocks the rewrite)*

Grep the bicomplex campaign configs for an explicit Hamiltonian threshold.
Whichever way it resolves, one of `research.tex:474-476` or
`postload_gate.py:22-23` is wrong and must be corrected.

### 4. PG-2 — make the boson-shell static/dynamic switch coherent *(S2)* — **DONE**

Option (b), enforced at the decode site where it cannot be bypassed: a candidate
that supplies no velocity key at all is static, full stop, whatever its
`shell_static` bit says. `resolved_shell_currents()` is now the single source of
truth for that decision and the decoder reads it. The space builder drops
`grtresna_shell_static` from any space that ends up with no velocity dimension —
note *ends up*: the v1 boson-shell builder strips the currents into `skip` and
then re-adds boson-specific ones when `static=False`, so the check has to run on
the finished dimension list, not on `skip`.

Recording: `unrecorded_shell_currents()` reports any non-zero current that
neither a search dimension nor an override declares, and the QD driver writes it
into the trajectory record as `unsearched_currents`. After the fix the stripped
case resolves to zero, so there is nothing left to hide — the field exists to
catch the next such default rather than this one.

### 5. PG-5 / PG-6 / PG-8 — mechanical *(S3)* — **DONE**

* PG-5: rather than derive the seed count from one expression, the worker RNG
  now mints one child stream per worker on demand from a live `SeedSequence`, so
  there is no count to get wrong. Successive `spawn(1)` calls return exactly the
  children `spawn(n)` would have, so existing workers' streams are unchanged.
* PG-6: `trajectory_lump(\d+)_` was added as a second pattern, but adding it to
  `_LUMP_KEY_RE` alone would **not** have been enough — the template then writes
  `grtresna_lump{k}_{suffix}` keys that a trajectory space does not contain, so
  it would still have been a no-op. The trajectory branch imposes the exotic
  flag only; the orbit is a searched dimension and the template has no business
  overwriting it.
* PG-8: `_bin_index` maps NaN to bin 0 and ±inf to the ends; the refill that
  keeps the pipeline full moved into a `finally`, so no exception from recording
  a result can cost an in-flight slot.

### 6. PG-7 — hygiene *(S4)* — **DONE**

Implemented as "search the half that survives decoding" rather than a magnitude
plus a separate sign dimension: when `trajectory_retrograde_only` is set, each
`trajectory_lump{k}_omega_rot` axis is narrowed to its non-positive half before
the search starts. Same effect, no new dimension, and no change to the stored
genome layout. The fold stays in the decoder as a backstop for pinned and legacy
values.

### 7. PG-9 — make the lump-count clamp loud, and reconcile the two caps *(S1)* — **1 and 2 DONE**

Three separate changes, in order of how much they buy:

1. **Refuse, don't clamp.** `SimulationParameters.hpp:126-127` should abort with
   the requested and available counts rather than silently truncating. A run
   that cannot do what its config asks should not start. Same for
   `rl_num_lumps` (`:66-67`) and `apply_rl_actions` (`RLActionApplier.hpp:31`),
   which clamp the same way.
2. **Let the wrapper know the cap.** The search space builder should read
   `GRTRESNA_MAX_INDEPENDENT_SCALARS` (or be told it) and refuse to declare
   trajectory dimensions it cannot drive, instead of emitting 56 and hoping.
3. **Decide whether 5 is still the right number.** The cap is a compile-time
   `static constexpr int` in `Source/Matter/GRTresnaScalarLayout.hpp:7` that also
   sizes the state vector. If 8-lump trajectories are wanted, this is the knob —
   raising it costs 6 more state components per extra lump-pair.

Re-run the 8-lump campaign's conclusions after (2); until then treat its archive
as a 35-dimensional search.

---

# PART C — CHECKED AND CLEAN (genome/decode segment)

* **Archive binning** — `_bin_index` clamps correctly for finite values; 1.0 →
  `bins−1`; verified 200/200 records in range, cells within [0,7], no duplicate
  eval ids, no NaN descriptors.
* **Pins wiring** — `--pin-dimension` validated against the space, removed from
  searched dims, merged into the base; per-eval `metadata.json` stores the full
  decoded overrides including pins (`core/evaluation.py:112-115`), so HQ replay
  physics matches the scored eval.
* **HQ replay path** — reads eval-dir metadata, deliberately bypasses
  `_vector_to_overrides` (`replay_eval.py:403-408`); replayed physics faithful.
  **This is why the 20.21% headline survives PG-1.**
* **QD↔CMA-ES encode ordering** — `_overrides_to_vector` is key-based, not
  positional; no off-by-one between drivers. Both build dims through the same
  `build_search_space`.
* **grtresna/GRTeclyn key partition** — `grtresna_*` filtered identically in all
  three eval paths; `trajectory_*` correctly flows to both.
* **Boson decode completeness** — every searched boson dimension is consumed
  (mode, layout, jitter, profile, exotic fraction/phase, axis, dipole,
  quadrupole, static, `grtresna_shift_seed`). No dead dimensions in
  shell/sh/ring/trajectory expansions **other than** PG-2's stripped velocities.
* **Exotic wedge selection** — ring inline copy and `_shell_exotic_ids` are the
  same algorithm.
* **Trajectory t=0 boson boost** — velocity is the exact spiral derivative,
  rotated by the same tilt matrix as position, capped at `boost_v_max`.
* **Resume counter / partial-dir scan** — `max(metadata, trajectory)`, no id
  reuse; prune protects in-flight and FTL-champion dirs.
* **Sampling weight normalization** — correct partial-sum draw.

---

# PART D — AUDITS STILL RUNNING

Five segments were dispatched in parallel; results land in this file as they
complete. Two have since been partly closed by the validators in PART F: the
`params.txt` key/arity sweep (row 1) and the GRTeclyn load side's component
contract (row 3) are now covered by executable checks — PG-9 came out of the
first, and the positional-vs-nominal gridinit split out of the second.

| segment | hunting |
|---|---|
| `GRTresnaConfig` → `params.txt` → GRTresna C++ | the bug-#19 shape (multi-token key read as scalar) on **every** vector key; dead keys; the μ=0 silent-default shape across mass/λ/μ₆/ω/κ |
| Chombo HDF5 → `.gridinit` | whether AMR levels are dropped the way the metric cache dropped them; χ = ψ⁻⁴ exponent; channel order; cell- vs node-centering |
| GRTeclyn load side | `ExternalGridInitialData` (two divergent copies); evolution potential vs solve potential term-by-term; trilinear-interpolation kink on refined levels |
| Acceptance gates | quantitative selectivity of PG-3; coverage map of gate × campaign; which drivers bypass the profile-contract rail |
| GRTresna solver internals | signed-lump source assembly; the exotic K=0 auto-switch self-consistency; convergence-norm dilution inside the solver; BC falloff vs GRTeclyn's Sommerfeld assumption |

---

# PART E — INTERACTION WITH `research.tex`

PG-1 adds a **new** required correction beyond the fifteen already listed in
`Debug.md`:

| tex lines | statement | status |
|---|---|---|
| 811-814 | "CMA-ES bought that by collapsing onto a quasi-static scaffold… The run is converged rather than merely lucky" | **WITHDRAW.** The near-zero rotation is the re-encode contraction, not a discovered optimum. Implied genome coordinate is 2.5% of the axis |
| 1099-1103 | "all five trajectories are retrograde but extremely slow… roughly an order of magnitude quieter than champion 87" | **WITHDRAW the interpretation.** The ~13.5× is the predicted `R0/v_max` contraction. The *values* are what was evolved and stay |
| 42-43, 1864-1867 | "QD search yields the champion" | **RESCOPE.** Valid for the 23 heritable dimensions; the 16 trajectory velocity/drift dimensions were random-sampled, not searched |
| 806-817 | CMA-ES basin / eight best genomes within 4 score points | **RESCOPE.** Consistent with a real basin *and* with ten dimensions frozen near zero. Cannot distinguish from stored data |
| 45-47, 1263-1283 | 20.21% headline advantage | **SURVIVES.** HQ replay does not re-convert; measured on the geometry `params.txt` specifies |

PG-9 adds one more, for any passage describing the 8-lump campaign:

| statement | status |
|---|---|
| the 8-lump trajectory search explored 8 lumps' trajectories | **RESCOPE.** GRTeclyn drives 5; lumps 5–7's 21 dimensions were inert. All 8 lumps exist in the initial data — only the trajectory driving is capped |

**The general rule for the rewrite:** PG-1 and PG-9 break *provenance* claims
(how the configuration was found), not *measurement* claims (what the
configuration does). Sort every sentence into one bucket or the other — the same
split that `Debug.md` §15 applied to the metric-cache bug.

---

# PART F — TESTS AND VALIDATORS ADDED

Per `Debug.md` §17's seventh habit — *"'physics unchanged' is not 'defect
removed' — every fix needs a test that fails before it and passes after"* — each
finding above now has an executable check. Run everything with:

```bash
cd grteclyn-wrapper && PYTHONPATH=src python -m pytest \
    tests/core/test_param_contract.py tests/core/test_cross_code_contract.py \
    tests/search/test_genome_roundtrip.py tests/search/test_shell_static_coherence.py \
    tests/search/test_pipeline_hygiene.py tests/search/test_retrograde_axis.py \
    tests/projection/test_postload_gate_composite_norms.py -q
```

### 1. Genome round-trip — `tests/search/test_genome_roundtrip.py` (20 cases)

**All green.** These are the PG-1 regression tests; each one failed before the
fix and passes after.

Asserts `E(D(x)) == P(x)` where `P = E∘D` is the projection the decoder actually
applies, that `P` is idempotent, that `P(x) == x` for every genome no cap binds,
and — for the one genome a cap does bind — the exact clamped value, so the
weaker-looking invariant is in fact the stricter test. Also: a noiseless mutation
of an elite reproduces that elite exactly from its stored genome, and a
strongly-rotating parent's children retain the signal rather than collapsing to
the 0.239 noise floor. A final test generalises to *every* searched dimension, so
a future transform that writes back into a genome key trips it without anyone
having to think of it.

*Note on the original plan's assertion.* `E(D(x)) == x` was specified in PART B
and is **not achievable**, because three of the decoder's transforms are genuine
physical projections rather than encodings. Relaxing to the projection form was
the correction, not a concession.

### 1b. The other five fixes

| file | covers | what would fail without the fix |
|---|---|---|
| `tests/search/test_shell_static_coherence.py` (10) | PG-2 | old decoder hands every lump `velocity = (0, 0, 0.35)` for a space that searches no velocity — asserted end-to-end on the lump list, not just the knob |
| `tests/search/test_pipeline_hygiene.py` (8) | PG-5, PG-6, PG-8 | `_bin_index(nan)` raised `ValueError`; the exotic template returned its input unchanged for trajectory dimensions; the pipeline refill sat after the recording instead of in a `finally` |
| `tests/search/test_retrograde_axis.py` (5) | PG-7 | the omega axis spans both signs while decode folds the positive half onto the negative |
| `tests/projection/test_postload_gate_composite_norms.py` (5) | PG-3 | the composite columns are not parsed at all; an all-zero refined-region norm from a single-level run reads as a clean pass |
| `tests/core/test_cross_code_contract.py` (27) | PG-10 | the C++ loader maps `.gridinit` columns by position while the Python side maps them by name, and nothing compares the two alias tables |

All seven files together: **102 passed**. Full wrapper suite: **915 passed, 6
failed, 10 skipped** (the 897 recorded here first, plus the 18 in
`tests/core/test_plotfile_scratch.py` — see below). Every failure is pre-existing
and unrelated — `h5py` and `cma` are not installed in this environment, and
`test_iterate` (2), `test_horizon_finder_guard` and
`test_qd_search_resume_continues_eval_counter` fail identically with these
changes stashed.

### 5. Node-local plotfile scratch — `tests/core/test_plotfile_scratch.py` (18)

Not a finding from this audit; recorded here because it changes where every
campaign writes and therefore what all of the above runs against. Plotfiles are
write-once/read-once/delete transients (~3.2 GB each at HQ resolution) and were
going to NFS for QD, CMA-ES, HQ and the post-load gate — which is what held the
node to two concurrent runs. `core/scratch.py` now maps every episode to a
node-local directory in `core/params.py::episode_path_overrides`, so no launcher
has to remember it. The promotion queues had done this by hand since
`b97868de`; this is the same policy at the one place every campaign passes
through.

The tests exist because moving the files breaks three things silently:
the scoring plotfile lookup (an episode that looks empty scores as though it
measured nothing), disk reclamation (deleting an unextracted plotfile from local
scratch loses the sample for good — on NFS it would have survived), and the test
suite itself (which was writing into the scratch root a live campaign uses).
Purging is ledger-gated on `small_data/consume_state.json`, and anything kept is
reported rather than silently retained.

**Verified on hardware, not only in tests.** A real GPU episode wrote 11
plotfiles, all to node-local disk, all extracted, with `small_data/`, frames and
`score.json` on the shared filesystem and **zero** heavy files left there. The
live `qball_traj_spin_v1` campaign shows the same at 12 evals.

### 6. Re-run of the PG-1 heritability check — IN FLIGHT

PART B item 1 ends: *"re-run one qball QD campaign and confirm the heritability
statistic stays near 0.5 rather than collapsing to 0.24."* That run is
`runs/grtresna_qd/qball_traj_spin_v1` (39-D, 200 evals, 8 GPUs, launched
2026-07-29), the first campaign in which the ten spin and drift axes are
heritable at all. Its dry-run configuration was checked before launch: 39
dimensions, all five `omega_rot` axes narrowed to `[-1, 0]` (PG-7), the raw
39-component genome stored in every trajectory record (PG-1), and no
`unsearched_currents` (PG-2). **The statistic itself is not yet measured** — use
the verification block in PG-1 once the run has ≳75 evals.

*Why the existing suite missed PG-1:* `test_trajectory_speed_cap.py` has 18
tests, all passing, all verifying the **forward** transform. Not one tests
whether it can be undone. The §17 "verified for the wrong property" trap.

### 2. params.txt ↔ C++ readers — `core/param_contract.py`, 28 tests

Static cross-check of an emitted `params.txt` against the ParmParse read sites
in both C++ trees. Catches four shapes, each pinned by a synthetic reproduction
of the bug that actually shipped:

| severity | check | bug shape |
|---|---|---|
| CRITICAL | vector emitted, every reader targets a scalar variable | #19 `recipe_scalar_field_signs` |
| CRITICAL | indexed key past `GRTRESNA_MAX_INDEPENDENT_SCALARS` | **PG-9** |
| CRITICAL | required physics constant not emitted, C++ silent default governs | #14 `recipe_scalar_mu = 0` |
| MAJOR | key emitted, read by neither C++ nor Python | dead knob |
| MINOR | destination type unresolvable | *checker gap, not an accusation* |

Two design points worth keeping:

- **Arity comes from the destination variable's type, not the function name.**
  `GRParmParse::load` is overloaded (`Source/IO/GRParmParse.hpp:27-72`): it
  forwards to the array `get`/`getarr` for `std::array`/`std::vector` and to the
  scalar `get` otherwise. A first version keyed off the function name and
  produced 8 false CRITICALs on correct code. Bug #19 was itself a *type* error
  — the sign array was declared `std::string` — so the type is the only honest
  discriminator.
- **Unresolvable is reported as MINOR, never CRITICAL.** A checker that cries
  wolf gets muted, and then it catches nothing at all.

Findings on the real tree dropped 37 → 3 as the scan learned about helper-built
keys (`load_coeff_array(pp, "recipe_chi_coeff", …)` builds `_0.._N`),
constructor-style declarations, and Python-side consumers. The 3 that remain:
`trajectory_r_min` and `trajectory_retrograde_only` are consumed by the Python
clamps (correct, now recognised); **`hdf5_subpath` is genuinely dead** — written
by `scripts/wormhole/run/wormhole_case.py:251`, read by nothing in either C++
tree or the wrapper.

### 3. GRTresna ↔ GRTeclyn ↔ gridinit — `core/cross_code.py`, 27 tests

The check the user asked for: *the values that go to GRTresna vs the values that
go to GRTeclyn*. Every run is configured twice, in two files, under two
different naming conventions, and nothing compared them until now.

| GRTresna | GRTeclyn | relation |
|---|---|---|
| `scalar_mass` | `recipe_scalar_mass` | equal |
| `scalar_lambda` | `recipe_scalar_lambda` | equal |
| `scalar_mu` | `recipe_scalar_mu` | equal |
| `lump<k>_exotic` | `recipe_scalar_field_signs[k]` | `sign = 1 − 2·exotic` |

A disagreement here does not raise: GRTresna converges, GRTeclyn runs, and the
initial data solves a *different Lagrangian* than the one evolving it — a
constraint violation present at t=0 by construction and indistinguishable from
physics. All 13 solved-and-evolved run pairs in `runs/` currently agree.

The `.gridinit` handoff is checked too, and it is the more dangerous of the two.
This is **PG-10**, and it is now fixed:

> **GRTeclyn's loader never read `component_names`.**
> `ExternalGridInitialData::compute` (`:69-91`) copied component `c` of the file
> into state variable `c`, for `c < min(num_components, NUM_VARS)`. The contract
> was purely **positional**. The Python analysis code indexes the same file **by
> name** (`comp_names.index("lapse")` and friends). Two codes, two different
> notions of what identifies a component, no cross-check. Reorder or insert one
> state variable and every field lands one slot off, silently — and
> `min(…, NUM_VARS)` truncated the tail without a word.

**Fix.** `Examples/RadialRecipe/ExternalGridInitialData.hpp` now resolves every
file column to a state-variable slot through `StateVariables::names`, iterates
destination slots rather than source columns, and refuses to load a file it
cannot map exactly — unknown name, duplicate name, or a `num_components` that
disagrees with the number of `component_names` (that last one decides the
payload stride, so guessing it misaligns every row after the first). A file
predating `component_names` still loads positionally but prints a warning
saying so. `Examples/RotatingWormholeCollapse/`'s copy already matched by name
and gained only the stride check.

**Verified as a pure safety net.** A standalone harness ran both mappings over
six real `.gridinit` files: `mapped=33 reordered=0` on every one. By-name and
by-position agree exactly on all existing data, so no run changes — the defect
is removed without the physics moving, which is the shape `Debug.md` §17 asks
for.

The validator parses the authoritative slot order out of `CCZ4StateVariables.hpp`
and `StateVariables.hpp` rather than duplicating it, and holds every real
`.gridinit` header to it. It also holds `grtresna/fields/boson_star.py:18-29` —
a hand-copied C++ enum whose only enforcement is the comment *"must match
StateVariables.hpp"* — to the actual enum.

Two documented aliases are accepted rather than flagged: `phi2`/`Pi2` and
`phi_lump0`/`Pi_lump0` are the same slots (`StateVariables.hpp:12-13` — the
single-complex-scalar model reuses the first lump slots). Flagging them would
have been wrong on every file in the tree.

### 4. Guard-the-guard tests

Most assertions above have the form "no finding of kind X". If a parser silently
regressed to returning nothing, all of them would go green while checking
nothing. Both suites therefore assert their inputs are non-empty: `>20` keys
parsed from a real `params.txt`, `>50` C++ read sites found, `>30` state
variables parsed, `names[0] == "chi"`.

### Still uncovered

- The domain mapping (`grtresna/domain/mapping.py`) — GRTresna's `L`/`N`/origin
  vs GRTeclyn's box and the gridinit's `dx`/`origin`. A half-cell offset here
  would be invisible to every check above. **This is the largest remaining gap.**
- The trilinear resample in `ExternalGridInitialData::compute` on refined levels
  (related to `Debug.md` #9).
- Whether GRTresna's boundary falloff matches GRTeclyn's Sommerfeld assumption.
