# Self-Gravitating Boson Star — Handoff (plain English)

## The short version

We were trying to drop a **self-gravitating boson star** into the simulation as a seed and have
it sit there as a stable, confined blob instead of dispersing. It kept spreading out. We tracked it
down to **four separate bugs** — one in the file loader, two in the physics of how the star was
painted onto the grid, and one in the "pump" that is supposed to hold the star in place. After
fixing all four and choosing a star on the *stable* branch, a single star now stays mostly
confined through a full evolution (at coarse resolution). It is **not perfect** — see the honest
caveats at the bottom — but it is dramatically better than before.

Nothing is committed yet. The work spans two repos: **GRTeclyn** (on branch
`feature/interstellar`) and the sibling **GRTresna** tree (`../GRTresna`). Both C++ trees were
rebuilt.

---

## What was actually broken (in plain English)

1. **The star's profile file was never being read.** The code that generates the star solves a set
   of equations and writes the answer to a text file (`qball_profile.dat`). The loader that reads
   that file back in choked on the `#` comment lines at the top, gave up immediately, and read
   **zero data points**. When it got nothing, it quietly fell back to a crude guess (a generic
   `sech` bump). So no matter how carefully we computed the real star, the simulation was always
   seeded with the wrong shape. This was the biggest problem, and it also explains why earlier
   "fixes" to the profile math seemed to change nothing.

2. **The star was painted as if it were sitting still in flat space.** A boson star is not static —
   its internal field rotates in a way that depends on how much the clock slows down near the star
   (the "lapse"). The code assumed the clock ran at normal speed everywhere (lapse = 1). Because of
   that, the momentum of the field was set wrong, which made the computed mass/energy wrong, which
   made the solved-for gravity wrong — so the star started out already out of balance and drifted.

3. **The "pump" was fighting the star instead of holding it.** There is a control mechanism (a
   feedback "trap") that nudges the field toward a target shape to keep it confined. The target was
   accidentally set to a value **10× too small** (it used a generic "well depth" number instead of
   the star's real central height). So the trap was constantly trying to squash the star down to
   a tenth of its size, which kicked off a wobble that eventually made it collapse.

4. **We picked a star that was doomed to collapse anyway.** Boson stars have a stable family and an
   unstable family, split by a peak mass (the "Kaup limit"). The value we were using put the star
   on the **unstable** side, meaning it collapses into a black hole on its own regardless of
   anything else. We moved to a value safely on the stable side.

---

## What we changed to fix each one

1. **Loader now reads the file properly.** It skips `#`/blank lines and accepts either 2 columns
   (flat-space case) or 3 columns (self-gravitating case). *(GRTresna: `BosonStarParams.hpp`.)*

2. **Star painted as a proper rotating star in curved space.**
   - The profile solver was rewritten to work in the coordinate system the evolution actually uses
     ("isotropic" coordinates) and now also outputs the lapse `α(r)` and the metric factor `ψ(r)`.
     *(`profiles/boson_star_ode.py`.)*
   - The file writer now emits a third column, the lapse `α(r)`. *(`solver/params.py`.)*
   - The painter now sets the field's momentum using the real, position-dependent lapse and uses
     `α(r)` as the starting lapse instead of 1. *(GRTresna: `ComplexScalarField.cpp`,
     `WriteOutput.H`.)*
   - We also corrected a factor-of-two in the energy formula so the Python solver and the C++
     evolution agree on the physics. This moved our computed Kaup limit from a wrong `0.435` to
     `≈0.62`, matching the textbook value `0.633`.

3. **Pump now aims at the right target.** We added a proper knob (`rl_pump_target_amp`) so the trap
   aims at the star's actual central height instead of the tiny "well depth" number, and the trap's
   width is now fitted to the star's real width. *(GRTeclyn: `RLMatterPumpParams.hpp`,
   `ComplexScalarField.impl.hpp`, `BiComplexScalarField.impl.hpp`,
   `RadialRecipeMatterDispatch.hpp`, `SimulationParameters.hpp`; wiring in
   `grtresna/matter/wiring.py`.)*

4. **Use a stable star.** For mass `m=1`, keep the central amplitude below the Kaup peak
   (`phi_c < 0.08`). We used `phi_c = 0.07`.

---

## What was tested

- **Unit / regression tests (Python):** full `tests/grtresna` suite passes (258 tests), and the
  code is lint-clean (`ruff`). New tests cover the isotropic solver, the lapse handling, and the
  pump-target wiring.
- **Solver sanity checks:** the profile solver reproduces the known Kaup limit (`M·m ≈ 0.62`) and
  the metric factors settle to their correct far-away values (`ψ→1`, `α→1`).
- **GPU smoke test (single star at rest, grid 128³, box size 64, run to t=16, coarse
  `max_level=1`):** with the stable star + fixed pump, the fraction of the field that stays
  confined climbs to ~0.97 around the middle of the run and is still ~0.90 at the end — compared to
  the old broken seed, which fell to ~0.58 and visibly dispersed. The run exits cleanly.

| time | old broken seed | fixed (stable star + fixed pump) |
|------|-----------------|----------------------------------|
| 0    | 0.96 | 0.74 (broad seed settling in) |
| 6.4  | 0.74 | **0.97** |
| 9.6  | 0.64 | **0.97** |
| 12.8 | 0.61 | **0.97** |
| 16.0 | **0.58** (dispersed) | **0.90** |

---

## Honest caveats (important — read before trusting the result)

- **The confinement number oversells it.** The "confined fraction" stays ~0.90 at the end, but that
  is measured against a fairly generous radius. A more sensitive measure, the field's RMS radius,
  tells a less rosy story: the star tightens to RMS ≈ 2.1 by t≈9.6 and then **spreads back out to
  RMS ≈ 4.2 by t=16** — i.e. it roughly doubles in size over the last third of the run. Looking at
  the raw movie frames, you can see this late-time spreading directly. So the star is **mostly held
  but slowly leaking / breathing outward**, not perfectly locked. The confinement metric is not an
  ideal detector of this; treat RMS radius (and the actual frames) as the better check.
- **High resolution still blows up.** At finer resolution (`max_level=3`) the run develops NaNs
  around t≈6–9. This is a numerical-relativity stability issue (constraint growth under
  regridding), separate from the seed/pump physics. The clean result above is only at coarse
  `max_level=1`.

---

## Suggested next steps

1. **Chase the late-time spreading.** Confirm from the frames + RMS radius how much the star is
   leaking, and decide whether it matters for the campaign. If it does, the most likely culprit is
   the pump *shape*: it is still a generic `sech` with the right height and width, not the star's
   exact profile. Loading the real tabulated `phi0(r)` as the pump target on the GPU would give an
   exact match and should reduce the drift.
2. **Fix the high-resolution instability** if long, fine runs are needed: add/raise Kreiss–Oliger
   dissipation, check the gauge (1+log lapse) start, tighten the CCZ4 damping (`kappa1/2`), and/or
   use fixed (non-regridding) refinement around the star.
3. **Pick a sensible default `phi_c`.** The current config default (0.08) sits right on the Kaup
   peak (marginal). Recommend a stable default (~0.06–0.07) for self-grav campaigns.
4. **Two-star case.** The per-lump plumbing (`total_alpha`, `total_pi2`, per-lump target amplitude)
   is already in place for a slow two-star orbit; not yet exercised.
5. **Commit both repos** once the above is settled (GRTeclyn `feature/interstellar` + GRTresna).

---

## How to reproduce

```
# Python side (grteclyn-wrapper/)
source .venv/bin/activate
python -m pytest tests/grtresna -q
ruff check src/grteclyn_wrapper/grtresna/profiles/boson_star_ode.py \
           src/grteclyn_wrapper/grtresna/matter/wiring.py \
           src/grteclyn_wrapper/grtresna/solver/params.py

# Rebuild C++ (both required)
#   GRTresna (initial data):
GRTRESNA_ENV=~/.mlspace/envs/grtresna CHOMBO_HOME=<...>/Chombo/lib \
  make -C ../GRTresna/Examples/BosonStarBH all -j8 MPI=TRUE DIM=3 OPT=HIGH DEBUG=FALSE
#   GRTeclyn (pump target; nvcc must be on PATH):
PATH=/usr/local/cuda/bin:$PATH make -C Examples/RadialRecipe \
  COMP=gnu USE_CUDA=TRUE USE_MPI=FALSE CUDA_ARCH=90 -j"$(nproc)"

# Winning single-star smoke test:
python grteclyn-wrapper/scripts/campaigns/hq/replay_eval.py runs/selfgrav_smoke/src_single_star \
  --name selfgrav_ml1_phic07_pumpfix --runs-dir runs/selfgrav_smoke/out \
  --gpu 0 --n-full 128 --l-full 64 --stop-time 16 --max-level 1 \
  --extra-override grtresna_bs_phi_c=0.07
# Check: small_data/confinement.dat  (col 5 = confined_frac, col 4 = rms_radius)
# and the movie frames under frames/phi_z/frames/
```

---

## Files touched

**GRTeclyn (`feature/interstellar`):**
`Source/GRTeclynCore/RL/RLMatterPumpParams.hpp`, `Source/Matter/ComplexScalarField.impl.hpp`,
`Source/Matter/BiComplexScalarField.impl.hpp`, `Examples/RadialRecipe/RadialRecipeMatterDispatch.hpp`,
`Examples/RadialRecipe/SimulationParameters.hpp`,
`grteclyn-wrapper/src/grteclyn_wrapper/grtresna/profiles/boson_star_ode.py` (new),
`grteclyn-wrapper/src/grteclyn_wrapper/grtresna/matter/wiring.py`,
`grteclyn-wrapper/src/grteclyn_wrapper/grtresna/solver/params.py`,
plus new tests `tests/grtresna/test_boson_star_selfgrav_ode.py`,
`tests/grtresna/test_selfgrav_profile_wiring.py`.

**GRTresna (`../GRTresna`):**
`Source/Matter/BosonStarParams.hpp`, `Source/Matter/ComplexScalarField.cpp`,
`Source/Tools/WriteOutput.H`.
