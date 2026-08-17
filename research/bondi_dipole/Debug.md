# Bondi dipole — review fixes & rerun plan

Updated 2026-08-17, after two external reviews of `bondi_dipole.tex`.
Constraint reminders: **this machine has no GPU** — every rerun below happens on
gigaclust2, and code reaches the cluster only via **commit → push → pull**
(never scp).

---

## 1. Answered: why the max level only triggers at t ≈ 48

Not a schedule — a consequence of the tagging criterion. The run executable
(`Examples/RadialRecipe`) tags with `ChiTagger`
(`Source/Tagging/ChiTagger.hpp`):

```
criterion = dx * sqrt( sum_ij (d_i d_j chi)^2 )  >= regrid_threshold (0.02)
```

i.e. it tags on the **curvature of χ**, not on the matter. The stars are
weak-field objects: |χ − 1| ≲ 0.02 spread over several m⁻¹, so at t = 0 the
criterion sits an order of magnitude below 0.02. `max_level = 1` was allowed
from the start, but nothing crossed threshold until the deepening canonical
well (χ_min 0.98 → 0.44 as the pair nears contact and radiation falls back)
steepened past it — at t ≈ 47.5 in the mixed runs, never in the singles.

Consequence: the entire science window (t ≤ 30) ran **unigrid at Δx = 0.5**,
~10 grid points per stellar rms radius. Both reviews flag this as the top
blocking issue. Fix is item 3-D below (field-based tagging from t = 0).

The paper now states this explicitly (Sec. Setup + Appendix).

## 2. Text fixes applied to bondi_dipole.tex (2026-08-17, no rerun needed)

- **Eq. (3)**: added the missing factor of 2 (`□Φ± = 2 ∂V/∂|Φ|² Φ±`), with a
  clause tying it to the ½-normalized kinetic term; `dV/dφ₀ ≡ 2φ₀ ∂V/∂|Φ|²`
  now defined next to Eq. (5).
- **Eq. (6) sign**: the flipped element was the Π *definition*, not the seed.
  The code convention is `∂ₜΦ = αΠ` (`complex_scalar.py: phase_velocity_pi2`,
  seeds `pi2 = −ω φ₀/α`), so the text now reads `Π = α⁻¹∂ₜΦ` and Eq. (6)
  keeps its minus sign — all three (ansatz, definition, seed) now consistent.
- **New theory paragraph** (Sec. II): the action
  S = ∫√−g[R/16π + ℒ(Φ₊) − ℒ(Φ₋)], Bianchi consistency, inherited strong
  hyperbolicity, and the named quantum-ghost/vacuum-decay objection
  (Carroll–Hoffman–Trodden 2003; Cline–Jeon–Moore 2004).
- **Forward link** added at "held together by something other than gravity …
  its own gravity pushes it apart" (the self-repulsion Forward identified).
- **"First" claims scoped**: "asymptotically flat, constraint-solved,
  dynamically evolved"; static phantom stars (Dzhunushaliev et al. 2014) and
  de Sitter negative-mass bubbles (Mbarek–Paranjape 2014) acknowledged;
  positive-mass-theorem citations (Schoen–Yau, Witten) added to the intro.
- **GPU-first claim** now shared with our own wormhole paper
  (arXiv:2604.00071), cited as companion work `\cite{shirokov2026}`.
- **Conclusion (3) demoted** to a theoretical expectation with the explicit
  P_signed(t) integral stated; abstract softened to match ("should carry …;
  direct momentum-budget measurement deferred").
- **Headline numbers**: 50–150× tied to the single-star noise floor with the
  ≥20× same-sign figure alongside (abstract, results, conclusions); 0.6% →
  "0.6–3.4% across the clean window"; 0.2% → "0.2–3%".
- **Table I** d₀ = 16 row: the wrong-signed canonical −0.34 now printed with a
  footnote explaining the boundary-loss bias.
- **Barycentre defined as an equation** (new Eq., Sec. Setup): field-weighted
  centroid with w_s = (|Φ_s|² + |Π_s|²)^{1/2} over the full domain
  (= accounting region), noted as coordinate/gauge-dependent.
- **Ψ₄ caveat**: recorded but not used quantitatively — both shells (R = 8, 16)
  sit inside the matter.
- **Superposition/phase paragraph** (Sec. III): lump phases 0 and π documented
  (inert for the mixed pair, physical for same-sector pairs); Helfer et al.
  CQG 39, 074001 (2022) cited; phase-offset PP rerun flagged as future work.
- **PM-eq defense**: why matching *ADM* mass (not bare mass/charge) is the
  right choice for Bondi's point-mass limit.
- **"Pull exceeds M/r²"** reworded to the point-mass value M/d² at the
  barycentric separation.
- **Table V** χ row relabeled "extremal χ (min for +, max for −)".
- **Δt = 0.02Δx** acknowledged as deliberately conservative and untuned.
- **Tagging wording** corrected ("χ-gradient" → curvature-of-χ with the exact
  criterion) + the t ≈ 48 explanation in Setup and Appendix.
- **AMR/gauge citations**: CTTK method paper (Aurrekoetxea–Clough–Lim 2023),
  1+log (Bona et al. 1995), Gamma-driver (Campanelli et al. 2006; Baker et al.
  2006), González–Guzmán–Sarbach 2009 (prior nonlinear ghost-scalar NR).
- **Title** trimmed: dropped "with GRTeclyn" (code stays in abstract/intro).
- **Draft box** expanded to 6 items and now points at this file.

Compile verified: 3 passes, zero errors, zero undefined citations/references.

## 3. Blocking items — need GPU reruns on gigaclust2 (priority order)

Every item here follows commit → push → pull; launch commands and env knobs
are in `results/bondi-dipole-runaway/LAUNCH.md` (single rank only — MPI
segfaults on the node).

- **A. Convergence series** (blocks *every* headline number): rerun `PM` and
  one control at N = 192 and 256. The configuration is symmetric under y→−y
  and z→−z, so quadrant symmetry buys ~4× and funds N = 256 at roughly
  current cost (~2–3 h/cell at N = 128 now). Richardson-extrapolate drifts,
  quote convergence order + error bars for 50–150×, 0.6–3.4%, ~10%.
- **B. Momentum-balance diagnostic** (blocks re-promoting conclusion 3): the
  scrutiny stream already exists — `GRTECLYN_SECTOR_DYNAMICS=1` (core
  positions, momentum balance, gauge check; see
  `grteclyn-wrapper/.../core/plot_consumer.py`). Run for `PM` and `PM-eq`;
  plot P_signed(t) = ∫(S_i[Φ₊] − S_i[Φ₋])dV decomposed by sector and by
  region (core spheres vs halo vs boundary flux). Must close to truncation
  error. Then restore conclusion (3) to a measurement.
- **C. Boundary problem** (blocks separation-scaling finals): implement a
  **sponge layer** (Kreiss–Oliger dissipation profile ramping up near the
  boundary) *and/or* enlarge the domain to L ≥ 128 (ideally 256 so boundaries
  stay causally disconnected past t = 60). Code change in RadialRecipe +
  **rebuild** (GPU build on cluster after pull). Move Ψ₄ shells out (R ≥ 32)
  or drop them. Massive-field modes reflect perfectly off Sommerfeld
  boundaries — this is why the window is capped at t = 60.
- **D. AMR from t = 0** (fixes §1): `ChiTagger` cannot see weak-field stars.
  Add a field-based tagger (tag on scalar amplitude/gradient of both sectors,
  or a fixed-region tagger around the cores) or drop the threshold to ~0.001
  (test for noise-tagging first). Target core resolution Δx ≤ 0.125 — two
  extra levels over the N = 128 base, or one over N = 256. Code change in
  `RadialRecipeLevel::tag_cells` + **rebuild**. Expect the ±8% breathing in
  Table V to shrink; re-quote it.
- **E. Core-weighted centroid rerun**: the diagnostic is "built and waiting"
  (paper, Next steps). Re-run the separation series (d₀ = 12, 16) with it;
  repair the d₀ = 16 canonical −0.34. Add a proper-distance separation
  between χ extrema along the axis as a gauge-robust check of Fig. 4c.
- **F. Noether-charge accounting**: per-sector conserved U(1) charge inside
  spheres — clean core-vs-halo bookkeeping the field-weighted rms stream
  can't give (it tracks the radiation bath).
- **G. Extra controls**: a phase-offset `PP` rerun (current pairs: phases 0/π)
  to bound the superposition systematic, and an **unequal-mass PP** control
  (mixed pair is unequal; an asymmetry-correlated diagnostic bias is
  currently undetectable).
- **H. Robustness point** at a second frequency or coupling (draft-box item).
- **I. CFL experiment** (cost saver, do first on one single-star cell):
  `dt_multiplier` 0.02 → 0.2 is ~10× cheaper; adopt only if the single-star
  health table reproduces.

Constraint targets after C+D: L∞(H) must stay ≪ 0.05 in the window (currently
0.05 at t ≤ 30, 0.9 by t = 60 — reviewers will reject on that curve; Fig. 8
needs to be re-made from the improved runs).

## 4. Local tasks (CPU, this machine — no GPU needed)

- **R(ω) check**: the paper claims higher ω → more compact stars; sextic
  solitonic stars can be non-monotonic in R(ω) (radius diverges in both the
  thick-wall ω→m and thin-wall ω→ω_min limits). Extract R(ω) from
  `analysis/star_family_scan.py` / `boson_star_ode`; if non-monotonic, fix the
  "more compact stars (higher ω)" sentence in Discussion.
- **Figures**: DONE 2026-08-17. `make_article_figures.py` restyled so line
  style, not gray shade, distinguishes runs: Fig. 4b — PM-eq now dash-dot
  (Φ₊) / dash-dot-dot (Φ₋) at 0.35 gray, controls dotted at 0.45 (was 0.78
  light gray); Fig. 4c — four distinct styles (solid / dash-dot / dashed /
  dotted); Fig. 5 — point-mass model dotted (canonical) / dash-dotted
  (phantom) at 0.4 (was 0.78); Fig. 8 — controls dotted / tightly dashed at
  0.45. Captions updated to match; regenerated and visually verified.
  The reviewer's "stray '1' artifacts in Figs. 1–8": **not reproducible** in
  the current PDFs (all eight inspected panel by panel) — likely fixed by the
  earlier "update figure annotations" commit or a viewer artifact on their
  side; re-check once more on the submission build.
- **Error bars**: derive preliminary scatter-based error bars from the four
  control runs for every headline number (superseded by A when it lands).

## 5. Bookkeeping for submission

- Remove the red draft box.
- Verify at proof stage: Hoffmann 1966 chapter title/page (1966 Hlavatý
  volume vs the 1965 King's College proceedings version); page numbers of
  Bonnor–Swaminarayan and Bičák et al. New bib entries added today were
  web-verified where possible (dzhunushaliev2014 = "Boson stars with
  nontrivial topology", arXiv:1409.6978); double-check gonzalez2009 II page
  015011 and mbarek2014 101502(R) formatting.
- "traversible" in the Shinkai–Hayward title is original — warn the
  copyeditor not to "fix" it.
- Sec. III D failure narrative stays in the main text; be prepared to move it
  to an appendix if a referee asks.
- After reruns A–E land: update Tables I/III/IV, Fig. 8, the clean-window
  ranges, and re-promote conclusion (3); only then consider the full-rewrite
  rule that applies to validated-fix campaigns.
