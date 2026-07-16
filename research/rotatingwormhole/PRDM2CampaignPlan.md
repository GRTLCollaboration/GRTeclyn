# PRD Natural-\(m=2\) Wormhole Campaign

The current matter setup is correct for an engineered, axisymmetric
support-removal experiment, but it is not yet a clean setup for a natural
gravitational-wave \(m=2\) instability. The PRD upgrade should first instrument
and validate the rotating Q-torus, then locate a genuinely unstable high-spin
branch, confirm seed-independent \(m=2\) growth, and finally repeat the
successful case with wave-zone extraction and rigorous controls.

## Scientific verdict

- The present `m=1` is the scalar phase winding in
  `Phi=f(rho,z) exp[i(m phi-omega t)]`; because `|Phi|^2=f^2` is axisymmetric,
  it does **not** imply GW `m=1` or `m=2`. The measured `>99.999%` GW power in
  `m=0` is therefore expected for the current collapse.
- The setup in [`article/research.tex`](article/research.tex) is valid as an
  engineered axisymmetric-collapse paper. It cannot support a claim that
  rotation causes the quadrupole.
- A natural bar mode cannot be “forced.” We can tune the system into an
  unstable regime and add only controlled numerical-noise probes. A genuine
  instability must show exponential matter-side `m=2` growth with a convergent
  rate independent of seed amplitude.
- Do **not** interpret a scalar `m_az=2` eigenstate or a quadrupolar kick as a
  bar mode: either can remain axisymmetric in density or impose the desired
  answer.

## 1. Correct and expose the actual matter model

- In
  [`../../Source/Matter/ComplexExoticScalarField.impl.hpp`](../../Source/Matter/ComplexExoticScalarField.impl.hpp)
  and
  [`../../Examples/RotatingWormholeCollapse/Main_SupportedWormhole.cpp`](../../Examples/RotatingWormholeCollapse/Main_SupportedWormhole.cpp),
  document mathematically that `support_strength` multiplies the phantom stress
  tensor while the Klein–Gordon RHS remains unchanged. State the resulting
  Bianchi/conservation violation explicitly.
- Resolve the manuscript/code mismatch: the implementation uses a linear ramp
  while [`article/research.tex`](article/research.tex) currently gives a cosine
  ramp. Prefer implementing the smooth cosine ramp and recording it in
  `params.txt`.
- Rename outputs and labels conceptually as `m_az` for scalar winding and
  `m_GW` or `m_density` for extracted modes.
- Audit AMR initialization: the current base-grid ID spacing is `0.5` while AMR
  evolves finer levels. For subtle nonaxisymmetric growth, generate or prolong
  initial data at adequate source-region resolution and measure level-by-level
  Hamiltonian/momentum residuals after initialization.

## 2. Add diagnostics before spending GPU time

- Extend
  [`../../Examples/RotatingWormholeCollapse/SupportedWormholeLevel.cpp`](../../Examples/RotatingWormholeCollapse/SupportedWormholeLevel.cpp)
  with normalized matter multipoles using a non-cancelling weight such as
  `|Phi|^2` or `|rho|`:
  - `C_m(t)=integral w exp(i m phi) dV / integral w dV`, for `m=1,2,3,4`.
  - Equatorial quadrupole `I_xx-I_yy`, `2I_xy`, bar amplitude, phase, and
    pattern speed.
  - Radial profiles of charge, angular momentum, and an effective
    angular-velocity/current profile for corotation tests.
- Add a postrun `bar_mode_analysis.py` that fits
  `|C_2|=A exp(gamma t)` only over an objectively selected linear-growth
  interval and reports `gamma`, pattern frequency, saturation amplitude, and
  uncertainties.
- Add matter-aware refinement or a fixed refined torus region alongside
  `ChiTagger`; first verify that changing tagging leaves the axisymmetric
  headline diagnostics unchanged within the existing resolution error.
- Validate the in-situ all-`m` Weyl extraction against the trusted Python path
  using an analytic/Teukolsky-wave test and the headline run. Dense in-situ
  output is needed because plotfile cadence is too sparse for growth rates and
  spectra.

## 3. Validate the Q-torus independently

- Keep `lambda=170`, `mu6=14450`, and `kappa=1` initially; changing the
  potential now would confound rotation, compactness, and support effects.
- Evolve each candidate Q-torus in 3D without a throat/gravity complication
  long enough to establish whether the tabulated eigenstate is stable and
  whether `C_2` remains at the resolution-dependent noise floor.
- Then run the throat configuration with constant support `S=1` to at least
  `t=80–100`. A viable bar-mode candidate must survive long enough for several
  e-foldings; the present rapid collapse after `t=8–10` gives noise-level
  asymmetries too little time to grow.
- Reject candidates whose apparent `m=2` amplitude tracks constraint growth,
  AMR regridding, boundaries, or seed amplitude linearly without exponential
  growth.

## 4. Search for a naturally unstable rotating branch

- Build constraint-solved ID families with the existing Q-torus solver,
  scanning **measured** quantities rather than assuming larger `omega` means
  larger spin. For Q-balls, `J=m_az Q`; `omega` also changes charge, radius,
  thickness, and compactness.
- Stage the search:
  1. `m_az=1`, three continuation points spanning low/mid/high compactness
     along the allowed `omega` branch.
  2. `m_az=2`, the same measured compactness range, as a higher-winding
     family—not as proof of a density bar.
  3. Only if needed, vary throat mass over a narrow bracket around the current
     value to move the torus-to-throat mass ratio.
- For every ID record `M_ADM`, `J_ADM`, `Q`, torus radius/thickness, peak
  `|Phi|`, constraint residuals, and a robust rotation measure. Do not use the
  fluid threshold `T/|W|=0.27` as universal for a phantom scalar; use it only
  as context.
- Pilot candidates cheaply at moderate resolution with constant support.
  Promote only cases showing at least three clean e-foldings of matter `m=2`,
  bounded constraints, and no simultaneous numerical pathology.

## 5. Prove that any \(m=2\) mode is natural

- For each promoted candidate run a seed hierarchy: exact symmetry/roundoff,
  `epsilon=10^-6`, and `epsilon=10^-5` density or field perturbations. The tiny
  seeds diagnose the mode; they must not determine its nonlinear amplitude.
- Required evidence:
  - The fitted growth rate `gamma` agrees across seed amplitudes and at three
    resolutions within uncertainty.
  - Onset times shift approximately as `-ln(epsilon)/gamma` while saturation
    amplitude and pattern speed agree.
  - `C_2` grows in the matter before the GW `Psi4^(2,2)` signal reaches the
    detector.
  - The `m=2` amplitude does not converge to zero; `m=1,3,4` are reported to
    exclude a different dominant instability.
  - The GW satisfies expected `+m/-m` relations and propagates consistently
    across extraction radii.
- Use two controls:
  - A physical `m_az=0` Q-ball ID matched as closely as possible in ADM
    mass/charge, rather than merely deleting the phase from the `m=1` torus.
  - A phase-stripped but externally supported toroidal-modulus control,
    explicitly labelled an engineered geometry control. Together they separate
    rotation from toroidal shape.
- After identifying a bar-unstable branch, test collapse in two ways:
  - Constant support followed by the same axisymmetric support withdrawal.
  - Preferably a constraint-solved `kappa<1` pump-free branch, which removes
    the time-dependent Bianchi violation from the final headline experiment.

## 6. Obtain PRD-level wave-zone evidence

- First run a dense finite-radius bridge on the proven `L=128`,
  `max_level=2`, MPI×2 setup with `R_ext={24,36,48}`, sponge `52–62`,
  sampling `Delta t<=0.25`, and stop time `>=80`. Require outer-radius overlap
  `>=0.98`, relative `L2` residual `<=15%`, and a stable three-radius `1/r`
  extrapolation.
- For the final bar-mode case, target `L=256`, `max_level=2`, MPI×4,
  extraction at approximately `R={60,80,100}`, sponge outside `R=104`, and
  evolution to `t≈150`. Confirm feasibility with a short initialization/memory
  test before committing the campaign.
- Wave-zone acceptance:
  - `|r Psi4^(2,2)|` varies by at most 5% across the outer radii after
    retarded-time alignment.
  - Outer-pair waveform residual is at most 5%, residual time shift at most
    0.1, and linear/quadratic `1/r` extrapolations agree within 5%.
  - The signal passes a three-resolution waveform comparison or a clearly
    separated Cauchy-error estimate.
  - Boundary/sponge contamination is causally excluded during the reported
    window.
- Ideal endpoint: add CCE and compare it with finite-radius extrapolation. If
  CCE is unavailable, keep the claim explicitly finite-radius and include
  extraction error bars; this is especially important for `m=0`, which is
  unusually gauge-sensitive.

## 7. Horizon, constraints, and reproducibility gates

- Before using “horizon formation,” add a genuine apparent-horizon finder and
  report horizon mass/spin and nonaxisymmetric horizon multipoles. Otherwise
  retain “trapped-surface proxy.”
- Report volume and throat-local `L2/Linf` Hamiltonian/momentum norms, CCZ4
  `Theta/Z_i`, and constraint convergence through bar growth and wave passage.
- Archive all promoted and control runs, exact ID files, dense mode data,
  parameter manifests, analysis scripts, code commit/tag, and
  figure-generation commands in a DOI-backed deposit.
- Expand the bibliography with spinning-boson-star bar-instability work
  (e.g. Di Giovanni et al., PRD 102, 124009), spinning Q-balls,
  Papaloizou–Pringle-type torus instabilities, and finite-radius versus CCE
  waveform systematics.

## 8. Manuscript decision gate

- If natural, convergent matter `m=2` growth is found, pivot the PRD paper to a
  nonaxisymmetric-instability result and use the present `m_GW=0` run as the
  low-spin/control case.
- If no candidate passes the growth tests, do not manufacture an `m=2`
  headline. Submit the narrower engineered support-removal paper after
  correcting the ramp equation, documenting controls, adding wave-zone
  evidence, and removing any implication that rotation sources the observed
  quadrupole.
- The publishable scientific outcome is valuable either way: a detected bar
  instability with a measured threshold, or a constrained null result over a
  documented rotating-Q-torus family.

## Execution checklist

- [ ] Align the support-ramp model, equations, terminology, and AMR
  initial-data treatment.
- [ ] Implement matter multipoles, bar growth/pattern diagnostics,
  matter-aware refinement, and validated dense all-`m` extraction.
- [ ] Validate candidate Q-tori in 3D isolation and under constant support
  before collapse runs.
- [ ] Construct and characterize `m_az=1` and `m_az=2` low/mid/high
  compactness ID families.
- [ ] Run a gated moderate-resolution search for exponential
  nonaxisymmetric growth.
- [ ] Confirm promoted `m=2` candidates with seed hierarchy, controls, and
  three-resolution growth rates.
- [ ] Run the dense `L=128` bridge and final `L=256` wave-zone campaign with
  finite-radius error estimates or CCE.
- [ ] Add horizon/constraint evidence, reproducibility archive, literature
  context, and revise the manuscript according to the outcome.
