# Capabilities of GRTeclyn

**GRTeclyn** is a modern C++ based code for doing numerical relativity simulations. It solves the Einstein equations of general relativity to evolve the metric components and matter fields from chosen initial conditions. It is built on [AMReX](https://amrex-codes.github.io/), a powerful framework for massively parallel block-structured adaptive mesh refinement (AMR), and supports both CPU and GPU architectures for high-performance computing.

GRTeclyn is the next-generation version of [GRChombo](https://github.com/GRTLCollaboration/GRChombo), with a modular and extensible design focused on clarity, maintainability, and performance portability. It is actively under development, and many of the core numerical relativity features from GRChombo are already implemented or underway.

Note that GRTeclyn is designed as a research tool, not a black box. Users are encouraged to modify and extend it for their physics cases. If you're a student or supervisor considering using GRTeclyn, please ensure your team has familiarity with numerical relativity, C++ and HPC environments.

---

## Current key features in the public code

### Formulation of the Einstein equations

* GRTeclyn uses the [CCZ4](https://inspirehep.net/literature/913488) formulation with optional constraint damping to evolve the Einstein equations.
* The implementation supports full 3+1D spacetimes with dynamic metrics and gauge evolution.
* The code is written modularly, with clearly separated evolution and diagnostic systems.

### Gauge choice and boundaries

* The moving puncture gauge is currently implemented by default, suitable for evolving black hole spacetimes.
* The gauge system is extensible and designed to support alternative slicing and shift conditions.
* The current version supports Sommerfeld (radiative) and reflective boundary conditions.

### Black hole spacetimes

* Binary black hole (BBH) simulations are supported using Bowen–York initial data.
* The **BinaryBH** example is fully functional. The code uses an approximate solution for zero spin binaries, or the TwoPunctures code for fully general initial momenta and spins.
* Puncture positions are tracked and recorded throughout the simulation using particle methods.
* Weyl extraction and modes are calculated using particle methods.

### Scalar fields on dynamical backgrounds

* **Oscillaton example** - for a minimally coupled scalar field, can easily be adapted to general potential.
* AMR interpolation tools — the scalar field example provides a simple example of how we enable variable interpolation and post-processing, based on particle methods.

### Adaptive Mesh Refinement (AMR)

* AMReX handles dynamic block-structured AMR with tagging criteria based on physical fields or custom diagnostics.
* Time sub-cycling and regridding are fully implemented and compatible with multi-level evolution.

### Finite differencing scheme

* The code supports 4th and 6th-order finite differencing in space and 4th-order Runge-Kutta in time.

### Checkpointing and restart

* Checkpointing is implemented for recovery after interruption.
* Simulations can be resumed from saved checkpoints with consistent memory and time-step restoration.

### Diagnostics

* Standard diagnostic variables (e.g., Hamiltonian and momentum constraints, stress-energy quantities) are implemented and can be output.
* A diagnostic system allows for easy addition of new variables to be monitored during the run. These are calculated only when required for output.
* Output is in AMReX-native plotfiles viewable in tools like **yt** or **AMReX-based viewers**, and small data files.

### Parallel performance and GPU support

* Designed for MPI-parallel runs using AMReX's distributed mesh management.
* GPU execution is supported using AMReX’s `GPU_LAUNCH` macros (e.g., CUDA/HIP/SYCL).
* OpenMP is supported for CPU parallelism on shared-memory systems.

---

## Features under development (coming soon)

These features are implemented in working branches or are actively being ported from GRChombo. They may be available on request (at your own risk pending testing).

* **Apparent horizon finder** — for locating black hole horizons in dynamical simulations.
* **Initial condition solver** — a link to the GRTeclyn solver GRTresna, to solve the Hamiltonian and Momentum constraints numerically for matter spacetimes.
* **ADM mass and momenta diagnostics** — for computing global quantities.
* **AMR reductions** - integration with AMReX methods for calculating integrals or max/min values across the adaptive mesh.

---

## Forthcoming features

The following features have been developed internally or in GRChombo and are expected to be ported or reimplemented in GRTeclyn. If you're interested in working on one of these, feel free to contact us — we may be able to share existing code for inspiration or adaptation.

* **Fixed background metrics** — support for scalar field evolution on Schwarzschild or Kerr backgrounds using analytic metrics - by adapting the existing [GRDzhaDzha](https://github.com/GRTLCollaboration/GRDzhadzha) code to GRTeclyn.
* **Cartoon formalism** — support for dimensional reduction using symmetries (e.g., axisymmetry).
* **Other matter models** — we will include perfect fluids and vector fields in future.
---

## Summary of development status

| Feature                                  | Status          |
|------------------------------------------|-----------------|
| CCZ4 evolution                           | ✅ Ported       |
| Adaptive Mesh Refinement (AMR)           | ✅ Ported       |
| Checkpointing & Restart                  | ✅ Ported       |
| Moving puncture gauge                    | ✅ Ported       |
| BinaryBH example                         | ✅ Ported       |
| TwoPunctures initial data                | ✅ Ported       |
| Diagnostic variables                     | ✅ Ported       |
| AMReX GPU offload                        | ✅ Ported       |
| MPI parallelism                          | ✅ Ported       |
| Scalar field matter model                | ✅ Ported       |
| Puncture tracking                        | ✅ Ported       |
| Particle-based diagnostics               | ✅ Ported       |
| Weyl scalar                              | ✅ Ported       |
| 6th-order finite differencing            | ✅ Ported       |
| Apparent horizon finder                  | 🔧 In progress  |
| Initial condition solver for matter      | 🔧 In progress  |
| ADM quantities                           | 🔮 Planned      |
| Cartoon formalism                        | 🔮 Planned      |
| CCE extraction                           | 🔮 Planned      |
| Fixed background metrics                 | 🔮 Planned      |
| Proca field matter and gauge fields      | 🔮 Planned      |
| Perfect fluid matter                     | 🔮 Planned      |

---
