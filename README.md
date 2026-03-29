# Wormhole Dynamics: Nonlinear Collapse and Gravitational-Wave Emission

This repository contains the numerical relativity C++ source code, simulation parameter files, extracted physical data, and Python analysis scripts for the paper:
**"Wormhole Dynamics: Nonlinear Collapse and Gravitational-Wave Emission"** by N. M. Shirokov.

This project is built on top of **[GRTeclyn](https://github.com/GRTLCollaboration/GRTeclyn)**, a GPU-accelerated numerical relativity framework (a port of GRChombo to the AMReX library).

## 🌌 About This Project

This repository provides the complete computational pipeline to simulate the 3D nonlinear dynamics of unstable Ellis–Bronnikov wormholes supported by a phantom scalar field. It includes the isotropic initial data setup, the forced quadrupolar collapse mechanism ($S_{\rm support}=0.5$, $A_\phi=+0.02$), and the extraction of the resulting gravitational-wave signals ($\Psi_4$).

### 📂 Navigation for Reviewers & Researchers

#### 1. C++ Simulation Code & Initial Data
The exact 3D Adaptive Mesh Refinement (AMR) setup, the initial data formulation, and the coupled Einstein–phantom-scalar evolution routines are located in:
* `Examples/SupportedWormholeCollapse/SupportedWormholeInitialData.hpp` (Isotropic Ellis-Bronnikov setup & perturbations)
* `Examples/SupportedWormholeCollapse/SupportedWormholeLevel.cpp` (CCZ4 evolution and diagnostic tracking)
* `Examples/SupportedWormholeCollapse/params_2gpu.txt` (The exact parameters used for the high-resolution production runs)

#### 2. Extracted Data (`.dat` files)
**Note on Data Storage:** To adhere to best practices and save space, the massive 3D `.hdf5` plotfiles (often hundreds of gigabytes) were discarded on the fly. All relevant physical observables were extracted during the evolution and are provided as lightweight 1D `.dat` files.
* Constraints, apparent horizon proxies, and collapse diagnostics can be found in `data_supported/data/` (or the equivalent extraction directories).
* The primary extracted gravitational-wave template is located at: `data_supported/small_data/psi4_mode_l2m0.dat`.

#### 3. Python Analysis & Visualization
The scripts used to generate the publication-ready figures (Richardson Convergence, Collapse Diagnostics, CWT Spectrograms) are modularized in the `src/visualisation/` directory.
* `src/visualisation/constraines/`: Plots the Hamiltonian and Momentum constraint convergence.
* `src/visualisation/extract_wave/`: Processes and plots the retarded-time $\Psi_4$ waveforms.
* `src/visualisation/make_evolution_panel/`: Generates the 2D $K_z$ trace snapshots of the phantom bounce.

#### 4. LIGO / GWOSC Data Search
If you are a gravitational-wave data analyst and wish to use this simulated wormhole collapse as a template in your own searches, we provide an automated search script:
* `src/search/search.py`
This script automatically downloads public LIGO data from the GWOSC API (e.g., GW190521), integrates $\Psi_4 \to h$, applies the correct physical scaling ($1000 M_\odot$ at $1$ Mpc), and executes a matched-filter search using `PyCBC`.

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