# Upstream merge plan — GRTLCollaboration/GRTeclyn:develop → `feature/interstellar`

**Written:** 2026-08-17
**Target branch:** `feature/interstellar`
**Author:** analysis session on the laptop — **nothing was compiled. Every claim below comes from source
inspection and `git merge-tree`, not from a build.** (AMReX does exist locally at
`/home/nik/Desktop/Wormhole/amrex`, so a slow local build is possible — see §10.)
**Validated:** 2026-08-17 — SHAs, divergence counts, conflict sets (both directions), commit split,
§6 quotes, §7 traps, §8 inventory, and §9 API re-checked against fresh fetches of both remotes.
**Purpose:** hand-off so the merge + port can be done on a machine with CPUs and AMReX.

---

## 1. TL;DR

* Syncing means taking **upstream PR #172**, which **replaced the entire matter interface**.
* Textual conflicts are small: **7 files** — identical set whether you merge into `develop` or
  `feature/interstellar` (verified both ways with `git merge-tree`).
* The real work is the **semantic port: 20 files, ~2,650 lines** of matter code that auto-merges cleanly
  and then **fails to compile**. Plus ~423 lines of now-obsolete `*Vars` helpers to delete, and 8 caller
  files to update.
* **Upstream ships no working example of the new matter interface** (§7.2). You are porting against an API
  with no reference caller. Biggest risk in the job.
* There is **no cheap middle path**: the 76 commits are one architectural change plus dependent work.

---

## 2. Reference SHAs

| Thing | SHA | Note |
|---|---|---|
| `upstream/develop` tip | `7dc7c8b5` | Merge PR #195 |
| `origin/feature/interstellar` | `48aa5a44` | **merge target** |
| `origin/develop` tip | `6c8e648a` | |
| merge base | `2bfd19ea` | common ancestor |
| **PR #172 merge commit** | `e5f8b380` | the architectural rewrite |

Divergence `feature/interstellar` ↔ `upstream/develop`: **499 ahead / 76 behind**.
(`feature/interstellar` is 73 ahead / 18 behind `origin/develop`.)

Commit split of the 76:

| Group | Count |
|---|---|
| Inside PR #172 (`Tensor`→`amrex::Array` + kernel fission) | **41** |
| After #172 (TwoPunctures #194, cleanup #204, lint #189, git-SHA #195, derivs fix #213, clang-tools #197) | **35** |

The 35 are **not** independent — e.g. `45b6996b BinaryBH: Update for tensor refactoring` assumes the new
API. Only a handful (git-SHA printing, lint config, clang-tools automation) are truly standalone.

---

## 3. State left on the laptop

* Added remote `upstream` → `https://github.com/GRTLCollaboration/GRTeclyn.git`
* Worktree `/home/nik/Desktop/github/GRTeclyn-merge`, branch `chore/merge-upstream`, based on
  **`origin/feature/interstellar`**, upstream-tracking unset so a stray `git push` can't hit a real branch.
* That worktree holds a **paused, conflicted merge**. It is disposable — recreate it, don't transfer it.
* `feature/interstellar` and the main working tree were **never modified**.

Discard laptop state:

```bash
git -C /home/nik/Desktop/github/GRTeclyn worktree remove ../GRTeclyn-merge --force
git -C /home/nik/Desktop/github/GRTeclyn branch -D chore/merge-upstream
```

---

## 4. Reproducing the merge on the build machine

Deterministic — same 7 conflicts every time.

```bash
cd <repo>
git remote add upstream https://github.com/GRTLCollaboration/GRTeclyn.git
git fetch upstream

git worktree add ../GRTeclyn-merge -b chore/merge-upstream origin/feature/interstellar
cd ../GRTeclyn-merge
git branch --unset-upstream          # avoid an accidental push to feature/interstellar

git merge upstream/develop           # → the 7 conflicts in §6
```

Escape hatches: `git merge --abort`, or `git worktree remove ../GRTeclyn-merge --force`.

Preview conflicts without touching any tree:

```bash
git merge-tree --write-tree --name-only origin/feature/interstellar upstream/develop
```

---

## 5. The architectural change (this is the whole job)

### Old (your code)

Driver pre-computes `Vars` / `D1Vars` / `D2Vars` / `AdvecVars` and passes them in:

```cpp
compute_emtensor(const Vars &vars, const D1Vars &d1,
                 const Tensor<2, amrex::Real> &h_UU,
                 const Tensor<3, amrex::Real> &chris_ULL) const;

add_matter_rhs(const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
               const D1Vars &d1, const D2Vars &d2, const AdvecVars &advec) const;
```

Plus **your** addition — a coords/time overload threaded through all three drivers:

```cpp
compute_emtensor(const Vars &vars, const D1Vars &d1,
                 const Tensor<2, amrex::Real> &h_UU,
                 const Tensor<3, amrex::Real> &chris_ULL,
                 const Coordinates &coords, amrex::Real time) const;
```

### New (upstream)

Matter class receives **raw grid indices + `Array4` + a derivative object** and computes its own
derivatives. `D1Vars`/`D2Vars`/`AdvecVars` are gone; `chris_ULL` dropped from `compute_emtensor`;
`Tensor<2,Real>` → `Tensor::Rank2`. The class is now templated on `deriv_t` as well.

```cpp
template <class potential_t = DefaultPotential,
          class deriv_t     = FourthOrderDerivatives>
class ScalarField
{
    using Vars = ScalarFieldVars;          // only Vars survives

    emtensor_t compute_emtensor(
        const int ix, const int iy, const int iz,
        const amrex::Array4<const amrex::Real> &state,
        const deriv_t &a_deriv,
        const Tensor::Rank2 &h_UU) const;

    void add_matter_rhs(
        const int ix, const int iy, const int iz,
        const amrex::Array4<amrex::Real> &rhs_state,
        const amrex::Array4<const amrex::Real> &state,
        const deriv_t &a_deriv) const;
};
```

Derivatives are fetched inside the body:

```cpp
const amrex::CellData<const amrex::Real> &state_cell_data = state.cellData(ix, iy, iz);
const Vars vars(state_cell_data);
auto d1_phi = a_deriv.d1_scalar(ix, iy, iz, state, c_phi);
```

### Why a "keep our version" shim is impossible

`CCZ4RHSWithMatter : public CCZ4RHS<gauge_t, deriv_t>`, and upstream replaced `CCZ4RHS`'s whole evaluation
model with index-based kernel fission (`compute_chi_and_h_ij(int ix,int iy,int iz,…)`,
`compute_A_ij_and_Theta_and_Gamma(…)`). Keeping your matter driver means rejecting #172 — the bulk of the
76 commits. Verified, not assumed.

---

## 6. The 7 conflicts, with resolution guidance

All 7 are the **same collision**: you thread `Coordinates &coords` + `m_time`; upstream removed that in
favour of indices + `Array4`. General rule: **take upstream's shape, then re-add `coords`/`time` as extra
trailing parameters.**

### 6.1 `Source/Matter/ScalarField.hpp`
Your side adds a coords/time `compute_emtensor` overload; upstream replaced the base signature.
→ Keep upstream's new signature, add your overload with `(…, const Coordinates &coords, amrex::Real time)`
appended. Keep `#include "Coordinates.hpp"`.

### 6.2 `Source/Matter/CCZ4RHSWithMatter.hpp`
Hunk 1: your `#include "Cell.hpp"` / `"Coordinates.hpp"` vs upstream removing them → keep `Coordinates.hpp`.
Hunk 2: `add_emtensor_rhs` signature, old Vars+coords vs new indices+Array4 → upstream shape + `coords`.

### 6.3 `Source/Matter/CCZ4RHSWithMatter.impl.hpp` — **most dangerous file**

Four hunks; two need real thought.

**(a) RHS zeroing — SILENT-WRONG-RESULTS TRAP.** Your side zeroes the RHS:
```cpp
for (int n = 0; n < rhs_cell_data.nComp(); ++n) { rhs_cell_data[n] = 0.0; }
```
Upstream's side does **not**, and states:
```cpp
// NB: the vacuum solution needs to be computed elsewhere!
// This will only compute the matter contribution
```
This is commit `0bfb08fa Matter: Split out vacuum solution from RHS`. **The matter driver no longer
computes the vacuum CCZ4 evolution.** Each `*Level.cpp` must now run the vacuum RHS and the matter
contribution and combine them. Take upstream's hunk without fixing every Level class and the code compiles
while silently producing wrong physics. **Highest-priority item in the port.** Composition recipe: §7.1.

**(b) your dispatch trait must survive:**
```cpp
if constexpr (std::is_same_v<matter_t, GRTresnaIndependentScalars> ||
              std::is_same_v<matter_t, ComplexScalarField> ||
              detail_matter::has_time_rhs_v<matter_t>)
{ m_matter.add_matter_rhs(rhs_cell_data, vars, d1, d2, advec, coords, m_time); }
else
{ m_matter.add_matter_rhs(rhs_cell_data, vars, d1, d2, advec); }
```
Re-express against the new signature:
`add_matter_rhs(ix, iy, iz, rhs_state, state, this->m_deriv[, coords, m_time])`.

### 6.4 `Source/Matter/ConstraintsWithMatter.hpp`
Include-only conflict. Keep `Coordinates.hpp`.

### 6.5 `Source/Matter/ConstraintsWithMatter.impl.hpp`
Upstream renamed the christoffel argument `d1` → `d1_h`:
```cpp
const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);
```
and calls `compute_emtensor(ix, iy, iz, state, m_deriv, h_UU)`. Take upstream, re-append `coords, m_time`.

### 6.6 `Source/Matter/Weyl4WithMatter.hpp`
Hunk 1: your `m_time` member init vs upstream's added
`// NOLINTEND(bugprone-easily-swappable-parameters)` → **keep both**.
Hunk 2: `add_matter_EB` signature → upstream shape (`ix,iy,iz,state`, `Tensor::Rank3`, `Tensor::Rank2`)
plus your `coords`.

### 6.7 `Source/Matter/Weyl4WithMatter.impl.hpp`
Matching call-site + definition changes for `add_matter_EB`. Same rule.

---

## 7. Traps

### 7.1 Vacuum/matter RHS split
See §6.3(a). Compiles fine, wrong answers. Affects every `*Level.cpp`. **Verify explicitly.**

The intended composition is demonstrated in `Tests/BSSNMatterTest/BSSNMatterTest.cpp` (~lines 138–170) —
this answers former open question §12.1. Per cell, in order:

1. `compute_chi_and_h_ij(ix, iy, iz, rhs, state)` — vacuum
2. `compute_A_ij_and_Theta_and_Gamma<formulation, covariantZ4>(…)` — vacuum
3. `calculate_gauge_rhs(…)` — vacuum
4. `CCZ4RHSWithMatter::operator()<formulation, covariantZ4>(…)` — matter contribution

Two landmines in that recipe:

* **Dissipation.** The vacuum path keeps `add_dissipation` in a separate `apply_dissipation` kernel
  (`Source/CCZ4/CCZ4RHS.impl.hpp`), which the matter composition **skips** — the matter `operator()` adds
  dissipation itself, once, at the end. A Level class that calls `apply_dissipation` *and* the matter
  operator applies dissipation twice; one that calls neither gets none.
* **Ordering.** Upstream's `add_emtensor_rhs` **assigns** (`=`, not `+=`) `rhs[c_Theta] = 0.0` in BSSN
  mode, so the matter kernel must run *after* the vacuum kernels, as in the test. (The assignment itself
  mirrors pre-split behaviour — BSSN doesn't evolve Theta — but it now stomps whatever the vacuum kernels
  wrote there.)

### 7.2 Upstream has NO working matter example — verify before trusting anything
* `Examples/ScalarField/` upstream still ships sources, but its only makefile is `GNUmakefile.old` →
  **not in the build**, and its `ScalarFieldLevel.cpp` still uses GRChombo idioms (`GRLevelData`,
  `BoxLoops::loop`, `MatterCCZ4RHS.hpp`). Dead code.
* `Examples/KleinGordon/` is a standalone wave solver — greps for `CCZ4RHS` / `add_matter_rhs` /
  `compute_emtensor` return **nothing**. Not a matter example.
* `Examples/BinaryBH/` is vacuum (+TwoPunctures).

The **only** references for the new interface are `Source/Matter/ScalarField.impl.hpp` and the tests
(`Tests/BSSNMatterTest`, `Tests/EMTensorTest`, `Tests/Weyl4WithMatterTest`). Read those first. The
vacuum + matter composition question is answered by `BSSNMatterTest` — see §7.1.

### 7.3 Dormant numerics bugfix
`57c4eb0e Derivs: Fix d1_tensor and d2_tensor` — `d1_tensor` wrongly assumed symmetric indices;
`d2_tensor` initialised `ivar{0}` instead of `ivar{ivar_0}`. **Dormant for you**: in your tree those names
occur only as local variables inside `Source/Grids/FourthOrderDerivatives.hpp` (6 hits, no external
callers). After the port, if matter classes start calling `a_deriv.d1_tensor(...)`, it goes live.

### 7.4 Revalidation
The port rewrites the numerics path of every matter class. Results must be re-validated against a
known-good checkpoint before any paper numbers are regenerated. **Do not treat the merged
`feature/interstellar` as paper-ready until that regression passes.**

---

## 8. Port inventory (from `feature/interstellar`)

Order: drivers → base `ScalarField` → derived matter classes → `*Vars` cleanup → Examples.
Line counts are `wc -l` on `origin/feature/interstellar` (an earlier draft measured the conflicted merge
worktree, which inflated the files touched by the merge).

### 8a. Drivers (do first; these are the conflicted files)
| Lines | File |
|---|---|
| 173 | `Source/Matter/CCZ4RHSWithMatter.impl.hpp` |
| 150 | `Source/Matter/ConstraintsWithMatter.impl.hpp` |
| 127 | `Source/Matter/Weyl4WithMatter.impl.hpp` |
| 92 | `Source/Matter/CCZ4RHSWithMatter.hpp` |
| 74 | `Source/Matter/ConstraintsWithMatter.hpp` (include-only conflict) |
| 67 | `Source/Matter/Weyl4WithMatter.hpp` |

### 8b. Matter classes — auto-merge clean, will NOT compile (**2,653 lines total**)
| Lines | File |
|---|---|
| 320 | `ControllerReservoirMatter.hpp` ← **only on `feature/interstellar`, absent from `develop`** |
| 203 | `ExoticScalarField.impl.hpp` |
| 178 | `BiComplexScalarField.impl.hpp` |
| 173 | `CCZ4RHSWithMatter.impl.hpp` |
| 158 | `ComplexScalarField.impl.hpp` |
| 150 | `ConstraintsWithMatter.impl.hpp` |
| 143 | `ComplexExoticScalarField.impl.hpp` |
| 135 | `GRTresnaIndependentScalars.impl.hpp` |
| 130 | `EMTensor.impl.hpp` |
| 128 | `GRTresnaIndependentScalars.hpp` |
| 127 | `Weyl4WithMatter.impl.hpp` |
| 122 | `ExoticScalarField.hpp` |
| 106 | `ScalarField.impl.hpp` |
| 99 | `EffectiveTeoMatter.hpp` |
| 92 | `BiComplexScalarField.hpp` |
| 91 | `ComplexScalarField.hpp` |
| 91 | `ComplexExoticScalarField.hpp` |
| 79 | `ScalarField.hpp` |
| 72 | `DustMatter.hpp` |
| 56 | `NoMatter.hpp` |

### 8c. `*Vars` helpers — likely **deletable** (423 lines)
`D1Vars` / `D2Vars` / `AdvecVars` no longer exist in the upstream interface:

| Lines | File |
|---|---|
| 66 | `ComplexScalarFieldAdvecVars.hpp` |
| 63 | `BiComplexScalarFieldAdvecVars.hpp` |
| 49 | `GRTresnaIndependentScalarsD1Vars.hpp` |
| 48 | `BiComplexScalarFieldD1Vars.hpp` |
| 46 | `GRTresnaIndependentScalarsAdvecVars.hpp` |
| 46 | `ComplexScalarFieldD1Vars.hpp` |
| 37 | `GRTresnaIndependentScalarsD2Vars.hpp` |
| 36 | `BiComplexScalarFieldD2Vars.hpp` |
| 32 | `ComplexScalarFieldD2Vars.hpp` |

Keep the plain `*Vars.hpp` (constructed from `CellData`) — those survive upstream.

### 8d. Callers to update (vacuum/matter split + signatures)
```
Examples/RadialRecipe/RadialRecipeLevel.cpp
Examples/RadialRecipe/RadialRecipeMatterConstraints.hpp
Examples/RadialRecipe/RadialRecipeMatterDispatch.hpp
Examples/RotatingWormholeCollapse/SupportedWormholeLevel.cpp
Examples/SupportedWormholeCollapse/SupportedWormholeLevel.cpp
Source/Grids/SpongeZone.hpp
Source/Matter/EMTensor.hpp
Source/Matter/MovingPunctureGaugeWithMatter.hpp
```

---

## 9. New derivative API (`Source/Grids/FourthOrderDerivatives.hpp`, upstream)

All take `(int ix, int iy, int iz, const amrex::Array4<const amrex::Real> &state, …)`.

| Method | Returns |
|---|---|
| `d1_scalar` | `Tensor::Rank1` |
| `d1_vector` | `Tensor::Rank2` |
| `d1_sym_tensor` | `Tensor::Sym12Rank3` |
| `d1_tensor` | `Tensor::Rank3` |
| `d2_scalar` | `Tensor::Sym12Rank2` |
| `d2_vector` | `Tensor::Sym23Rank3` |
| `d2_tensor` | `Tensor::Sym34Rank4` |
| `d2_sym_tensor` | `Tensor::Sym12Sym34Rank4` |
| `advection` (+ scalar/vector/tensor overloads) | scalar / `Tensor::Rank1` / `Tensor::Rank2` / `Sym12Rank2` |

`SixthOrderDerivatives` also exists upstream — matter classes are templated on `deriv_t`, so both must
instantiate.

---

## 10. Build

Possible on the laptop — AMReX sits at `/home/nik/Desktop/Wormhole/amrex`, with `g++ 13.3.0`, `mpicxx`,
and 8 cores — but slow; nothing was compiled during analysis. Preferred: the build machine.

```bash
export AMREX_HOME=/path/to/amrex     # on the laptop: /home/nik/Desktop/Wormhole/amrex
cd Tests && make -j$(nproc)          # port the tests FIRST — smallest surface, real reference
cd ../Examples/<target> && make -j$(nproc)
```

Suggested order:
1. Resolve the 7 conflicts (§6).
2. Build `Tests/` and fix until green. `BSSNMatterTest`, `EMTensorTest`, `Weyl4WithMatterTest` are your
   reference implementations of the new interface.
3. Port matter classes one at a time (§8b), rebuilding tests after each.
4. Port the Level classes (§8d) — mind §7.1.
5. Delete the obsolete `*Vars` helpers (§8c).
6. Build the Examples.
7. Regression-run against a known-good checkpoint before trusting any number.

---

## 11. Merging back

```bash
# only once tests are green AND a regression run matches a known-good checkpoint
git checkout feature/interstellar
git merge chore/merge-upstream
```

Merging `develop` is a separate, optional exercise — it produces the same 7 conflicts but does **not**
exercise `ControllerReservoirMatter.hpp`, so a green `develop` build would not prove the port complete.

Do **not** rebase — 499 commits replayed over a changed interface means resolving these conflicts hundreds
of times. Merge is the only sane option.

---

## 12. Open questions for upstream

1. ~~How are Level classes meant to compose vacuum RHS + matter contribution after `0bfb08fa`?~~
   **Answered** — `Tests/BSSNMatterTest` demonstrates the composition; recipe and landmines in §7.1.
2. Is `Examples/ScalarField` intentionally dead, or is a ported version coming?
3. Is there a migration guide for out-of-tree matter classes? You are almost certainly not the only
   downstream fork hit by #172.
