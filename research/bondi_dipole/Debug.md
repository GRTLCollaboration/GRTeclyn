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
- **C. Boundary problem + wave-zone Weyl** (blocks separation-scaling finals
  AND the radiated-energy question raised in the new Sec. VI): implement a
  **sponge layer** (Kreiss–Oliger dissipation profile ramping up near the
  boundary) *and/or* enlarge the domain to L ≥ 128 (ideally 256 so boundaries
  stay causally disconnected past t = 60). Code change in RadialRecipe +
  **rebuild** (GPU build on cluster after pull). Massive-field modes reflect
  perfectly off Sommerfeld boundaries — this is why the window is capped at
  t = 60. Weyl program on the enlarged domain (the near-zone analysis of
  2026-08-17 cannot settle any of these):
  - Move the Ψ₄ shells into the wave zone (R ≥ 32, several shells to verify
    the r⁰ peel of r·ψ₄; today's shells at 8/16 show r·ψ₄ falling 5–10×
    between them — pure near zone).
  - **Extend the extraction to l = 3 (and l = 4)** in
    `grteclyn-wrapper/.../extraction/psi4.py` — currently l = 2 only. The
    odd-l/even-l interference carries the momentum-flux beaming along the
    runaway (x) axis; l = 1 does NOT exist for ψ₄ (spin weight −2), which the
    paper now states explicitly against the "dipole GW" expectation. This is
    a **local CPU code change** that must land before the rerun (extraction
    happens live in the streaming consumer; plotfiles are purged).
  - **Net radiated energy sign** for the mixed pair (positive, negative, or
    inter-sector cancellation) from wave-zone dE/dt — the headline follow-up
    number; Sec. VI promises it.
  - **ADM mass + momentum surface integrals** at the outer boundary vs time
    (should stay ≈ const; PM-eq should give ADM ≈ 0 globally) — decouples the
    conservation statement from scalar-field centroids.
  - **Bonnor–Swaminarayan benchmark**: compute the B&S acceleration parameter
    for (M₊, M₋, d₀) and compare with the measured early-time coordinate
    acceleration, per run — the strong-field analogue of Table II.
  - Optional figure for the follow-up: shift vector β^x profile along the
    axis (frame dragging of the runaway — matter-free velocity proxy).
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
0.05 at t ≤ 30, 0.9 by t = 60 — reviewers will reject on that curve; the
constraints figure (`fig_constraints.pdf`, appendix — Fig. 9 since the Weyl
figure took Fig. 8) needs to be re-made from the improved runs).

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
- **Weyl analysis**: DONE 2026-08-17. New Sec. VI ("The gravitational-wave
  signature") + `fig_weyl.pdf` (new `fig_weyl()` in
  `make_article_figures.py`), from the packed `psi4_mode_l2_all.dat` streams
  (l = 2 spin-weighted modes, R = 8/16). Headline numbers (all reproducible
  by rerunning `fig_weyl()`, which prints them):
  - A_l2 at R = 16 grows 297× (PM), 832× (PM-eq), 481× (MP mirror) above the
    t ≤ 10 floor; PP bounded near 1e-3 with no secular ramp; MM irregular
    (its own spreading).
  - corr(A_l2, |Ẋ|) = 0.995 over 5 ≤ t ≤ 55 (pair-centre speed, ±2 central
    difference); fitted κ = 0.016.
  - Mirror reproduces the PM amplitude history to 4.7% (max-norm).
  - m = ±1 carry 0.6% of the l = 2 power (x-axisymmetric pattern about the
    z-decomposition axis).
  - Spectrum of A_l2 is secular: 99.7% below ω = 0.4; 4e-4 in the ω = 0.55
    band, 8e-6 at 2ω — NOT scalar-bath leakage.
  - Near-zone verdict: A(R=8)/A(R=16) = 5.0 at t = 30, 9.9 at t = 60 →
    r·ψ₄ ~ R⁻²…⁻³, nothing wave-zone; hence no energy/flux claims in the
    paper — that's rerun item C.
  - Physics note recorded in the paper (dipole question): signed source is
    conserved ⇒ signed-dipole current = conserved signed momentum ⇒ no
    mass-dipole radiation at leading order (point-mass reciprocity
    M₊a₊ = |M₋|a₋); ψ₄ (spin −2) has no l = 1 mode at all. Counterpoint to
    the naive "negative masses ⇒ dipole GW" reading of Trivedi–Loeb
    (arXiv:2605.10976), which itself notes universal coupling evades dipole
    bounds.
  - Raw waveform added as `fig_weyl` panel (a) (2026-08-17): monotone
    chirpless ramp in PM vs oscillatory wander in PP/MM at R = 16.
  - **WARNING — wormhole pipeline on this data**: running
    `process_wave/plot_extracted_psi4.py --combined` on
    `data/pair_pm/psi4_mode_l2m0.dat` works mechanically but most panels
    are physically invalid on near-zone data: the strain/LIGO panel reports
    "445 Hz, SNR 133 at 10 Mpc, horizon 166 Mpc" and E_rad = 1.9e-4 M —
    these assume wave-zone radiation and MUST NOT be quoted anywhere; the
    propagation-speed panel returns v = ∞ (no burst front to track — the
    secular ramp peaks at end-of-run at both shells simultaneously); the
    spectrogram is empty (all power below its f = 0.1 floor). Only the
    waveform and PSD panels are meaningful. After rerun item C (wave-zone
    shells), the full pipeline becomes valid — rerun it then, including the
    strain/LIGO comparison the wormhole paper made.

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
- Sec. III D failure narrative **removed from the article** (2026-08-17, user
  decision: debug text, not article text). Full content preserved in §6 below
  and in git history. Cross-references repaired: contribution (v) now points at
  Caveats only; the couplings, dressed-star, and seeding paragraphs keep the
  physics conclusions ("first couplings at which a dressed ω=0.55 star
  exists", "any rescaling de-tunes the eigenstate") without the campaign
  story. The two load-bearing residual facts survive elsewhere: the 30×
  tighter-solve probe lives in the Caveats resolution paragraph; the t = 0
  seed verification lives in contribution (iii).
- 2026-08-17: Sec. VI (GW signature) added with `fig_weyl.pdf`; abstract,
  contributions (iv), roadmap, conclusion (1), and Next steps updated to
  match; Sec. II strengthened (Gleiser–Dotti exponential modes + no-puncture
  argument; Bonnor's inverted fluid interior, verified against the PDF;
  companion wormhole e-folding 0.11 M). New bib entry `trivedi_loeb2026`
  (arXiv:2605.10976, web-verified 2026-08-17); Table V labels shortened to
  fix a 20.6 pt overfull. Figure numbering shifted: Weyl figure is Fig. 8,
  constraints figure is now Fig. 9.
- After reruns A–E land: update Tables I/III/IV, the constraints figure
  (`fig_constraints.pdf`), the clean-window ranges, and re-promote
  conclusion (3); only then consider the full-rewrite
  rule that applies to validated-fix campaigns.

## 6. Removed from the article: the Sec. III D failure narrative (for the record)

Cut 2026-08-17. Kept here because two of the four root causes are the reason
the t = 0 verification gate exists, and future reruns must not regress them.

Before any runaway could be claimed, a single lump had to hold its own size;
four campaigns failed that gate. Heavy lumps collapsed (min χ → 0.009 by
t ≈ 25), light lumps evaporated (rms 5 → 30 by t ≈ 50), and intermediate
lumps did both — indicting the *shape* of the seed, not its weight. Four root
causes, in order:

1. **The amplitude clamp.** The seeding script clamped the amplitude to the
   flat-space thin-wall value √(3λ/4μ), 4.5% below the solved eigenstate: an
   exact Q-ball shape at the wrong height is not a Q-ball — it breathed, shed,
   and evaporated. *Fix:* seed at the solved φ_c.
2. **No dressed equilibrium existed.** The corrected flat-space seed still
   dispersed: a flat-space eigenstate is an equilibrium only in flat space,
   and the dressed family showed the target star did not exist at the original
   couplings. *Fix:* weaken the couplings to the working values and seed the
   dressed profile with its own lapse (paper Eq. 6).
3. **Two seeding paths disagreed silently.** The constraint solve used the
   tabulated star while the evolution's matter reconstruction fell through to
   a Gaussian: the solved metric backed one object, the evolution started
   another. *Fix:* both paths resolve the same star, enforced to 1e−10 by a
   regression test.
4. **The phantom path had never run.** An old veto downgraded exotic lumps to
   canonical matter whenever self-gravity was on — correct for a
   gravity-*bound* star, irrelevant for one bound by self-interaction.
   *Fix:* lift the veto. The first phantom star evolved matched its predicted
   profile, with the predicted mirror-image metric signature (χ just *above*
   unity).

Residual effects attributed rather than fixed: the ±8% radial breathing is
resolution-intrinsic (re-solving with 30× tighter momentum residuals
reproduced it crest for crest); the rms-radius stream tracks the shed
radiation bath, not star health.

The exportable lesson: **verify the seed before believing the evolution.**
Every run was checked at t = 0 against the independently integrated star
solution; two of the four root causes were caught by that check, not by
watching movies.
