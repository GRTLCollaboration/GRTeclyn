# Wormhole Dynamics: Nonlinear Collapse and Gravitational-Wave Emission

This repository contains the numerical relativity C++ source code, simulation parameter files, extracted physical data, and Python analysis scripts for the paper:
**"Wormhole Dynamics: Nonlinear Collapse and Gravitational-Wave Emission"** by N. M. Shirokov.

Preprint available at [arXiv:2604.00071](https://arxiv.org/abs/2604.00071).

This project is built on top of **[GRTeclyn](https://github.com/GRTLCollaboration/GRTeclyn)**, a GPU-accelerated numerical relativity framework (a port of GRChombo to the AMReX library).

## About This Project

This repository provides the complete computational pipeline to simulate the 3D nonlinear dynamics of unstable Ellis–Bronnikov wormholes supported by a phantom scalar field. It includes the isotropic initial data setup, the forced quadrupolar collapse mechanism ($S_{\rm support}=0.5$, $A_\phi=+0.02$), and the extraction of the resulting gravitational-wave signals ($\Psi_4$).

### Project layout

#### Modified simulation (C++ / parameters)
The executable is built from **`Examples/SupportedWormholeCollapse/`** (`GNUmakefile`, `Make.package`). **Example-specific** sources include `SupportedWormholeInitialData.hpp`, `SupportedWormholeLevel.hpp` / `.cpp`, `Main_SupportedWormhole.cpp`, `SimulationParameters.hpp`, `StateVariables.hpp`, `PhantomDecayPotential.hpp`, and **`params_2gpu.txt`**.

**Shared matter code** is compiled from **`Source/Matter/`** (declared in `Source/Matter/Make.package` and pulled in via `src_dirs` in the example `GNUmakefile`). Files there:

| Role | Headers |
|------|---------|
| Scalar field (φ, derivatives, exotic branch) | `ScalarField.hpp`, `ScalarField.impl.hpp`, `ScalarFieldVars.hpp`, `ScalarFieldAdvecVars.hpp`, `ScalarFieldD1Vars.hpp`, `ScalarFieldD2Vars.hpp`, `ExoticScalarField.hpp`, `ExoticScalarField.impl.hpp` |
| CCZ4 + matter (RHS, constraints, stress tensor, Weyl) | `CCZ4RHSWithMatter.hpp`, `CCZ4RHSWithMatter.impl.hpp`, `ConstraintsWithMatter.hpp`, `ConstraintsWithMatter.impl.hpp`, `EMTensor.hpp`, `EMTensor.impl.hpp`, `Weyl4WithMatter.hpp`, `Weyl4WithMatter.impl.hpp` |
| Gauge & potential | `MovingPunctureGaugeWithMatter.hpp`, `DefaultPotential.hpp` |

The same `GNUmakefile` also pulls in other **`Source/`** trees (CCZ4, Wormholes, GRTeclynCore, …); initial data for *this* run is implemented in the example as **`SupportedWormholeInitialData.hpp`**, not `Source/Wormholes/WormholeInitialData.hpp` (that header serves other wormhole setups). The example evolution includes matter headers directly, e.g. `CCZ4RHSWithMatter.hpp`, `ConstraintsWithMatter.hpp`, `Weyl4WithMatter.hpp`, and `ExoticScalarField.hpp` from `SupportedWormholeLevel.cpp`, and `PhantomDecayPotential.hpp` uses `ScalarFieldVars.hpp` from **`Source/Matter/`**.

#### Shell automation — `grteclyn-wrapper/scripts/`
| Script | Role |
|--------|------|
| `plot_run.sh` | While the run writes plotfolders: `consume_plotfiles` **extracts** the useful observables ($\Psi_4$, optional frames) into **`.dat`** files and **PNG** frames; processed plotfile trees can be **deleted** so raw simulation dumps are not kept—**disk optimisation**. |
| `plot_diagnostic.sh` | After a run: constraint norms, collapse diagnostics, $\Psi_4$ panels into `grteclyn-wrapper/src/grteclyn_wrapper/visualisation/plots/`. |
| `move_files.sh` | Copy key run data + figures into **`grteclyn-wrapper/output/SimResults/<run_folder>/`**. |

Legacy shims under `src/scripts/` forward to the paths above.

#### LIGO / GWOSC matched-filter search — `grteclyn_wrapper.gw_search`
Template-based search in public strain data (GWOSC via GWpy, matched filtering with PyCBC). Entry point:

```bash
cd grteclyn-wrapper && uv run python -m grteclyn_wrapper.gw_search.main
```

Defaults use an extracted waveform under `grteclyn-wrapper/output/SimResults/…/psi4_mode_l2m0.dat`; override with `--data-path`. Methodology and options: **`grteclyn-wrapper/src/grteclyn_wrapper/gw_search/README.md`**.

#### Archived simulation outputs — `grteclyn-wrapper/output/SimResults/`
Packaged results per run (parameters, extracted `.dat` data, diagnostic plots, etc.), populated by `move_files.sh` and comparable archiving.

#### Plotting & analysis — `grteclyn_wrapper.visualisation`
Full usage, CLI flags, and workflows: **`grteclyn-wrapper/src/grteclyn_wrapper/visualisation/README.md`**. Subfolders at a glance:

| Path | Purpose |
|------|---------|
| `visualize/` | 2D field slices from plotfiles; optional MP4 (`ffmpeg`). |
| `visualize/make_movies.py` | Stitch existing PNG frame folders into MP4. |
| `make_evolution_panel/` | Multi-panel strip figure from saved frame PNGs. |
| `extract_wave/` | Extract $\Psi_4$ (mode $l{=}2$, $m{=}0$) from plotfiles; waveform + PSD. |
| `process_wave/` | **`consume_plotfiles`** (stream plotfiles → `psi4_mode_l2m0.dat`, optional frames) and **`plot_extracted_psi4`** (plots from `.dat`, optional strain / LIGO-style overlays). |
| `diagnostic/` | Collapse diagnostics from `collapse_diagnostics.dat` (optional `areal_radius.dat`). |
| `constraines/` | $L_2$ norms of Hamiltonian and momentum constraints. |
| `figures/` | Standalone schematic publication figures (not tied to one dump). |
| `plots/` | Default output location for `plot_diagnostic.sh`. |

#### Extracted data (`.dat` files)
**Note on data storage:** Large 3D `.hdf5` plotfiles are often discarded; observables are extracted during the run as lightweight 1D `.dat` files.
* Constraints, apparent horizon proxies, and collapse diagnostics: `data_supported/data/` (or the extraction directories for your run).
* Primary gravitational-wave template example: `data_supported/small_data/psi4_mode_l2m0.dat` (paths may differ if you use another output directory).

---

## 🛠️ About the Underlying Framework: GRTeclyn

GRTeclyn (previously referred to as GRAMReX) is a new numerical relativity code developed by the[GRTL Collaboration](https://www.grtlcollaboration.org) that is currently still under development. It is a port of the [GRChombo code](https://github.com/GRChombo/GRChombo) (based on the Chombo libraries) to the [AMReX](https://amrex-codes.github.io/) library in order to leverage AMReX's support for GPUs and ongoing active development.

The AMReX documentation can be found [here](https://amrex-codes.github.io/amrex/docs_html).

The name follows a similar pattern to GRChombo, namely "GR" + "\<Tool in a foreign language\>". In this case, "Teclyn" is a Welsh word for "Tool".

### Development status
Please consult this [documentation page](https://grtlcollaboration.github.io/GRTeclyn/#summary-of-features) for a list of the development status of specific features.

### Documentation
Documentation can be found [here](https://grtlcollaboration.github.io/GRTeclyn/) (under construction). Note that the GitHub wiki is no longer in use.

The documentation contains useful information on obtaining and building the code, prerequisites, and running the binary black hole example.

### License
GRTeclyn is licensed under the BSD 3-Clause License. Please see [LICENSE](LICENSE) for details.
