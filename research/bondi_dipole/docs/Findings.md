# Findings

Things measured during the campaign that are not run bookkeeping: artefacts
identified and named, effects explained, questions closed. The run-book and the
gates live in `GPU_RUN_PAPER.md`; this file holds the "why" answers — the ones a
referee asks after looking at a figure or a movie.

Each entry states what was seen, what it is, and whether it touches the result.

---

## The flat pulse in the shift movies is emitted by the outer boundary at t=0 — it is not a reflection

**Measured 2026-08-22** on `runaway_pair_d10_L64_N256_lev0` (L=64, N=256,
`max_level=0`), from the cached frame slices covering t = 0 → 16.8.

### What is seen

The `shift1` z-slice movie shows a faint straight vertical line near the right
edge of the box that walks steadily inward. It is easy to read as radiation from
the pair bouncing off the wall. It is not.

Tracking the line frame by frame along y = 32:

| t | x of the front |
|---|---|
| 1.6 | 63.12 |
| 4.8 | 58.12 |
| 8.0 | 53.62 |
| 12.8 | 46.88 |
| 16.8 | 41.12 |

A straight-line fit over t = 2.4 → 16.8 gives speed **1.389** moving in −x, and
back-extrapolates to **x = 64.45 at t = 0** — the right wall, to within two
cells. The left wall launches the mirror image: it is first resolved apart from
the star edge at t = 7.2 (x = 23.12) and reaches x = 14.12 by t = 13.6, speed
1.406, and the *inward* branch from x = 0 arrives at x ≈ 23.9 by t = 16.8,
matching the 1.39 prediction.

So both x-walls emit a planar front at t = 0 and it travels inward.

### Why it cannot be a reflection

Three independent reasons, any one of which is sufficient:

1. **Nothing has reached the wall yet.** The stars sit at x = 27 and x = 37, so
   the nearest wall is 27 away; the fastest mode in the system travels at
   ≈ 1.41. The earliest anything the pair emits can *touch* a wall is t ≈ 19,
   and the frames showing the front start at t = 0.8. There is nothing out there
   to reflect.
2. **It is superluminal.** 1.39 in code units, where c = 1. Physical radiation
   cannot do this; it is consistent with the 1+log lapse mode, whose speed is
   √(2/α) ≈ 1.41 where the lapse has relaxed to ≈ 1 in the far field. The 2%
   deficit is within the tracking error of locating a peak to a cell over a
   14.4-long baseline.
3. **It is flat, not curved.** The front is a plane spanning the full y range.
   Anything originating at the pair is a spherical arc, and those arcs are
   visible in the same frames as expanding circles around the stars.

### What it actually is

Start-up junk in the gauge sector. The initial data hands the evolution a shift
of exactly zero everywhere, which does not satisfy the condition the outer
boundary enforces. The boundary corrects itself on the first step and the
correction propagates inward at the gauge speed. This is the ordinary reason
the first light-crossing time of a run is not trusted; it is worth naming here
only because the movie makes it look like something else.

### Why only the shift movie shows it

It is in every field, at the same place, moving at the same speed. Measured from
the frame-to-frame change, y-averaged, at t = 8.0:

| field | front sits at | pulse amplitude | field's own scale | fraction of the picture |
|---|---|---|---|---|
| `shift1` | x = 53.4 | 1.9e-06 | 2.0e-05 | **10%** |
| `chi` (χ−1) | x = 53.4 | 3.1e-06 | 1.0e-02 | 0.03% |
| `local_speed` | x = 54.6 | 3.2e-06 | 1.0e+00 | 0.0003% |

The shift is a nearly empty field at this stage, so a 2e-06 pulse occupies a
tenth of its colour range. χ carries real structure a few hundred times larger,
so the identical pulse is a rounding error on its bar and never renders. The
contrast is a property of the colour scale, not of the physics — and it is why
"only one movie shows it" is not evidence that only one field is affected.

### Does it touch the result

It reaches the pair at t ≈ 23 and has passed beyond them by t ≈ 46. Both
crossings fall inside the completed base rung, so the question is answerable
from data already in hand.

`runaway_pair_d10_L64_N128_lev0`, Hamiltonian constraint norm:

| t | 9 | 23 | 45 | 69 |
|---|---|---|---|---|
| Ham | 4.03e-06 | 4.49e-06 | 5.87e-06 | 7.01e-06 |

Monotone and smooth — no step, no kink, no bump at either crossing time. The
slow rise across the whole window is the ordinary secular drift of the run, not
a response to the front. The fitted acceleration over the same window agrees
with the late-time value (1.428e-04 at t = 63 against 1.440e-04 over
t = 157 → 199), so the pass of the front leaves no measurable trace on the
quantity the paper reports.

**Conclusion: cosmetic. It is visible in one movie because that movie has an
empty colour range, and it is invisible in every diagnostic that matters.**

### How to reproduce

The cached frame slices are enough — no plotfiles, no re-run:

```
.venv/bin/python - <<'PY'
import numpy as np, glob
d = 'runs/bondi/staging/runaway_pair_d10_L64_N256_lev0/bondi_sg_pair_pm_eqm_w075'
xs = np.arange(256) * 0.25 + 0.125
for f in sorted(glob.glob(f'{d}/frames/_slice_cache/shift1_z/*.npz')):
    z = np.load(f, allow_pickle=True)
    a = z['arr'].astype(float)
    prof = np.abs(a - a.mean()).mean(axis=0)     # y-average: planes survive, arcs wash out
    R = xs > 40
    j = np.argmax(prof[R])
    print(f"{float(z['time']):6.2f}  x={xs[R][j]:6.2f}  amp={prof[R][j]:.2e}")
PY
```

Swap `shift1_z` for `chi_minus_1_z` or `local_speed_z` to see the same front at
the same x with an amplitude that never rises above their own noise.
