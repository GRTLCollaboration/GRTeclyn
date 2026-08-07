# How the matter was made to stop dispersing

Before any runaway could be claimed, a single lump had to hold its own size.
Four campaigns failed that gate. This is the diagnosis trail — useful for the
article's methods section, and as a catalogue of failure modes for anyone
seeding solitons into a constraint solver. The full chronological record is
[`debug_log/bondi_dipole_debug.md`](debug_log/bondi_dipole_debug.md).

---

## 0. The failures, and what their *pattern* said

| campaign | λ, μ | core amplitude | outcome |
|---|---|---|---|
| strong | 640, 85 333 | 0.075 | **collapse** — min χ → 0.009 by t ≈ 25–30 |
| weak | 2 560, 1 365 333 | 0.0375 | **dispersal** — rms 5 → ~30 by t ≈ 50 |
| weak v2 (corrected seeds) | 2 560, 1 365 333 | 0.0375 | dispersal on the *identical* curve |
| mid | 1 280, 341 333 | 0.053 | **both at once** — core crunch *and* envelope blow-off |

Heavy lumps collapsed, light lumps evaporated, and the geometric-mean rung did
both simultaneously. No amplitude rung could fix that, which said the problem
was not the *weight* of the blob but its *shape*: whatever the solver was
producing was far from any self-holding soliton at these couplings, so it
violently re-arranged regardless of mass.

Two early hypotheses were killed by measurement, which is what made the real
cause findable:

- **"The seed knobs are wrong."** Setting `bs_phi_c` and `bs_profile_width` to
  their textbook values changed the t = 0 sector integrals **not at all** —
  identical to 9 significant digits. Those knobs are only the solver's initial
  *guess*; it relaxes to its own state. Dead end, cheaply proven.
- **"The pair overlap tears the lumps apart."** A *lone* lump, with nothing to
  overlap, dispersed on exactly the pair's curve. Refuted.

## 1. Root cause A — the amplitude clamp (a 4.5 % error that mattered enormously)

Both painters (C++ and Python) rescale the tabulated Q-ball profile by
`amp / φ_c`. The campaign path set `amp` through `cap_well_depth`, which clamps
to the **thin-wall estimate** √(3λ/4μ) — but the true eigenstate amplitude φ_c
sits above it. At the campaign couplings the clamp landed at

```
amp / φ_c = 0.95452
```

so **every seed was born at 95.45 % of its own eigenstate amplitude** — an
exact Q-ball shape at the wrong height, which is not a Q-ball at all. It
breathed, shed, and (with gravity too weak to hold it) evaporated.

*Fix:* `grtresna_qball_exact_amplitude=1` — opt-in, so existing campaigns keep
bit-identical initial data — and, on the self-gravitating path, the writer
back-writes `amp := solved φ_c`, making the rescale exactly 1.

*Verification:* the fixed seed reproduced its predicted t = 0 fingerprint
(measured 36.2998 / 5.2209 against a predicted 36.300 / 5.221).

## 2. Root cause B — no dressed equilibrium existed at those couplings

The fix worked, and the lump still blew off. That falsification was the
important one: an exact flat-space eigenstate is only an equilibrium *in flat
space*. Once gravity is switched on, the object must sit on the
**gravitationally dressed** family — and solving that family at fixed frequency
showed the target simply was not on it:

> At the weak rung, the dressed family near ω = 0.55 maps to feather-light
> stars (ADM ≈ 0.05) or ultra-compact unstable ones (ADM ≈ 2.3) — nothing near
> the flat ball's E ≈ 0.28. **The lump had no equilibrium to settle into.**

Scanning rungs found the first one where the dressed ω = 0.55 star exists:
λ = 10240, μ = 21845333 ("ultraweak"), giving φ_c = 0.019695, ADM = +0.0640,
α(0) = 0.977. The family scan in `stars/star_family.csv` puts the branch edge
at ω ≈ 0.52–0.53 for this rung.

*Fix:* move to the ultraweak rung **and** seed the dressed star itself —
including its own lapse — rather than a flat-space profile.

*Verification:* the canonical star survived to t = 40 with an intact core; the
old seeds had unravelled by t ≈ 12.

### 2.1 A solver trap worth recording

Shooting on ω to find the star fails on this branch. Heavy-branch sextic stars
sit a fraction of a percent below the effective-potential top, so at fixed φ_c
the eigenvalue in ω is an exponentially thin needle: the shooter steps over it
and lands on the light branch. Observed concretely — asking for the weak-rung
ω = 0.55 star returned a **5× lighter ω = 0.75 star**. The working
parameterisation bisects φ_c at fixed integration-frame ω and outer-iterates ω
so the rescaled ω/α_∞ matches the request
([`MATTER_MODEL.md`](MATTER_MODEL.md) §3.2).

## 3. Root cause C — two painters that disagreed silently

With the dressed star seeded, a subtler inconsistency surfaced: the constraint
solve used the tabulated star, but the Python **repaint** that builds the
evolution's initial data fell through to a Gaussian envelope for the
self-gravitating profile type. The solved metric therefore backed one object
while the evolution started from another — a mismatch invisible in any single
file, and exactly the kind that shows up as "the blob relaxes at t = 0".

*Fix:* the repaint dispatches on the same profile type and resolves the same
star (same cache entry, same `gravity_sign`) as the solver.

*Lesson generalised:* when two code paths must paint the same object, test that
they agree **pointwise**, not that each is individually plausible. The
regression test asserts the repaint reproduces the star table to 1e-10 and is
*not* the Gaussian.

## 4. Root cause D — the phantom sector had never actually been run

The exotic (phantom) path carried a veto: with self-gravity on, an exotic lump
was silently downgraded to canonical matter, on the reasoning that a
gravity-*bound* phantom star cannot exist (its self-gravity is repulsive). True
for the mini-star path — but not for a soliton bound by its own
self-interaction, where gravity is only a dressing.

*Fix:* lift the veto when sextic couplings and a target frequency are present,
and solve the phantom star with `gravity_sign = −1` (metric sources only). Emit
its own table, since sharing the canonical one would seed it off-equilibrium.

*Verification:* the first lone phantom star ever solved **and** evolved here hit
its predicted fingerprint (20.989 / 5.427 against 20.99 / 5.43), showed the
predicted mirror-image metric signature (χ just **above** 1, where the canonical
star sits below), and turned out to be the most stable object in the campaign.

## 5. Two residual effects, correctly attributed

Both were candidate "bugs" that measurement reassigned:

**Breathing (±8 %) is resolution-intrinsic, not residual-driven.** The obvious
suspect was a dirty constraint solve — the first dressed run exited at Mom
residual 0.64 % (the gate was 1 %), and that residual radiated as a visible χ
ring that sloshed the envelope. Tightening the gate to 0.1 % / 0.002 gave
Mom = 0.018 % (30× cleaner) and killed the ring — **and the breathing was
unchanged, tracking the dirty run crest-for-crest.** It is discretisation of the
continuum star on N = 128 / L = 64. The lever is resolution, not tolerance.

**The rms stream is not a star-health diagnostic here.** It climbed 5 → 13 while
the core sat intact and breathing: it tracks the *shed radiation bath*
(amplitude-linear weight × r² leverage), not the star. Core health has to be
read from peak amplitude, confined fraction and min χ. This also explains the
"dispersal" language in the early campaigns — some of it was bath growth.

**The bath cannot leave the box.** Massive-field radiation against
massless-wave (Sommerfeld) boundaries reflects back and washes over the star;
canonical-sector χ slides monotonically (0.99 → 0.85 by t = 40) while the
phantom sector's does not. This is what caps mixed cells at t = 60 and is the
main obstacle to longer runs — a sponge layer is the fix.

## 6. Process lessons

1. **Verify the seed before believing the evolution.** Every cell in this pack
   was checked at t = 0 against an independently computed fingerprint. Two of
   the four root causes above were caught by that check, not by watching movies.
2. **Falsify cheaply, in the right order.** The two dead-end hypotheses cost one
   measurement each; had they been "fixed" speculatively, the amplitude clamp
   would have stayed hidden behind plausible-looking changes.
3. **A fix that does not work is information.** Exact amplitude *not* solving
   dispersal is what proved no dressed equilibrium existed — the second root
   cause was only reachable through the first fix's failure.
4. **Run controls under identical numerics.** The null results in `pair_pp` and
   `pair_mm`, and the single-star noise floor (drift ≈ 0.06), are what turn "the
   blobs moved" into "the mixed pair, and only the mixed pair, moved".
5. **Gate risky seams before they can contaminate a result.** The never-run
   phantom path was evolved *alone* first; had it been debuted inside the pair
   cell, a phantom-side bug would have been indistinguishable from the physics
   under test.
6. **Stop detached campaigns with the sanctioned tool.** Killing a worker
   advances the orchestrator to the next cell; the stop script kills the
   orchestrator first, sweeps workers, then verifies.
