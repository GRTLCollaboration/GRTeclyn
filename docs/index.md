# Welcome to GRTeclyn's documentation!

GRTeclyn is a flexible numerical relativity code based on [GRChombo](https://github.com/GRTLCollaboration/GRChombo) and
built on the [AMReX](https://amrex-codes.github.io/amrex/) library. The
GRTeclyn source code is available on
[GitHub](https://github.com/GRTLCollaboration/GRTeclyn) under the
BSD-3-Clause license.

If you can't find the answers to your questions here, we also have a GitHub [Wiki](https://github.com/GRTLCollaboration/GRTeclyn/wiki) and Slack channel (contact us for details on how to sign up). 

## Where should I start?

Here are some suggestions depending on your goal:

### *I want to know if GRTeclyn is the right code for my research project*

Start with the [**Capabilities**]() page for an overview of GRTeclyn's current features.

GRTeclyn is **not** a black-box code. You’ll need to understand how it works, and may need to write your own tools for specific applications. If you're a supervisor planning to assign a student to work on the code, please ensure you have sufficient background in **numerical relativity** and **C++/GPU computing** to support them.

For simpler projects e.g. for masters students, consider using the 1D code [**engrenage**](https://github.com/GRTLCollaboration/engrenage/wiki), which is a python based code that runs easily on laptops.


#### *I want to compile and run an example (e.g. binary black hole merger)*

Follow sequentially the **Getting started** section. This step-by-step guide assumes only basic command-line experience (that of most new PhD students). We also assume access to a high-performance computing cluster (preferably with GPU support). While we provide laptop tips, most realistic simulations require HPC.


#### *I already have NR experience and want to use or adapt the code*

1. Still begin with the **Getting started** section to ensure your environment is set up properly.
2. Then head to the section **Doing Physics with GRTeclyn**, which covers the code structure, design philosophy and guidance for adding your own physics

If you have questions, feel free to [**Contact us**](Contact-us) — we’ll try to help. But do keep in mind our unofficial motto:

> _“If you want it, you build it.”_

Before publishing results, don’t forget to check the [**License**](License) and [**Citation**](Citation) pages.


#### *I found a bug*

Please file an [issue](https://github.com/GRTLCollaboration/GRTeclyn/issues) on **GitHub**:

- Include your code version/commit hash
- Try to include as much information as you can, including: 
  - Detailed steps to reproduce the issue
  - Your build/runtime environment
  - Any error messages
  - The expected outcome/output of your run

You can also [**Contact us**](Contact-us), but GitHub Issues are preferred for code-related problems.


#### *I want to contribute to GRTeclyn*

That's great! See [**Contributing to GRTeclyn**](Contributing-to-GRTeclyn) for contribution guidelines and tips to get started.


###️ Disclaimer

> **GRTeclyn is a research code under active development.**  
> We aim to provide a stable core, but not all features are guaranteed to work perfectly.  
> **Use it at your own risk!**

## Capabilities of GRTeclyn

**GRTeclyn** is a modern C++ based code for doing numerical relativity simulations. It solves the Einstein equations of general relativity to evolve the metric components and matter fields from chosen initial conditions. It is built on [AMReX](https://amrex-codes.github.io/), a powerful framework for massively parallel block-structured adaptive mesh refinement (AMR), and supports both CPU and GPU architectures for high-performance computing.

GRTeclyn is the next-generation version of [GRChombo](https://github.com/GRTLCollaboration/GRChombo), with a modular and extensible design focused on clarity, maintainability, and performance portability. It is actively under development, and many of the core numerical relativity features from GRChombo are already implemented or underway.

---

### 🚀 Current key features in the public code

#### Formulation of the Einstein Equations

* GRTeclyn uses the [CCZ4](https://inspirehep.net/literature/913488) formulation with constraint damping to evolve the Einstein equations. 
* The implementation supports full 3+1D spacetimes with dynamic metrics and gauge evolution.
* The code is written modularly, with clearly separated evolution and diagnostic systems.

#### Gauge Choice and boundaries

* The moving puncture gauge is currently implemented by default, suitable for evolving black hole spacetimes.
* The gauge system is extensible and designed to support alternative slicing and shift conditions.
* The current version supports Sommerfeld (radiative) and reflective boundary conditions.

#### Black Hole Spacetimes

* Binary black hole (BBH) simulations are supported using Bowen–York initial data. 
* The **BinaryBH** example is fully functional for boosted cases. Currently the code uses an approximate solution but we plan to integrate the TwoPunctures code shortly, which will allow fully general initial momenta and spins.
* Puncture positions are tracked and recorded throughout the simulation using particle methods.

#### Scalar Fields on dynamical backgrounds

* Oscillaton example under development but minimally coupled scalar field is working.

#### Adaptive Mesh Refinement (AMR)

* AMReX handles dynamic block-structured AMR with tagging criteria based on physical fields or custom diagnostics.
* Time sub-cycling and regridding are fully implemented and compatible with multi-level evolution.

#### Finite Differencing Scheme

* The code supports 4th-order finite differencing in space and 4th-order Runge-Kutta in time.
* Support for higher-order stencils (6th order) is planned in future releases.

#### Checkpointing and Restart

* Checkpointing is implemented for recovery after interruption.
* Simulations can be resumed from saved checkpoints with consistent memory and time-step restoration.

#### Diagnostics

* Standard diagnostic variables (e.g., Hamiltonian and momentum constraints, stress-energy quantities) are implemented and can be output.
* A diagnostic system allows for easy addition of new variables to be monitored during the run. These are calculated only when required for output.
* Output is in AMReX-native plotfiles viewable in tools like **yt** or **AMReX-based viewers**, and small data files.

#### Parallel Performance and GPU Support

* Designed for MPI-parallel runs using AMReX's distributed mesh management.
* GPU execution is supported using AMReX’s `GPU_LAUNCH` macros (e.g., CUDA/HIP/SYCL).
* OpenMP is supported for CPU parallelism on shared-memory systems.

---

#### 🛠️ Features under development (coming soon)

These features are implemented in working branches or are actively being ported from GRChombo. They may be available on request (at your own risk pending testing).

* **Gravitational wave extraction** — Weyl scalar calculations for waveform analysis, and output of data for CCE extraction.
* **Apparent horizon finder** — for locating black hole horizons in dynamical simulations.
* **Initial condition solver** — a link to the GRTelcyn solver, to solve the Hamiltonian and Momentum constraints numerically for scalar or vacuum spacetimes. Integration of TwoPunctures as discussed above.
* **AMR interpolation tools** — to enable conservative variable interpolation and post-processing, based on particle methods.
* **ADM mass and momenta diagnostics** — for computing global quantities.
* **Higher-order spatial stencils** — planned support for 6th-order differencing.

---

### 🔮 Forthcoming features

The following features have been developed internally or in GRChombo and are expected to be ported or reimplemented in GRTeclyn. If you're interested in working on one of these, feel free to [[Contact us | Contact-us]] — we may be able to share existing code for inspiration or adaptation.

* **Fixed background metrics** — support for scalar field evolution on Schwarzschild or Kerr backgrounds using analytic metrics - by adapting the existing [GRDzhaDzha](https://github.com/GRTLCollaboration/GRDzhadzha) code to GRTeclyn.
* **Massive vector (Proca) fields** — a matter class for evolving massive vector fields.
* **Cartoon formalism** — support for dimensional reduction using symmetries (e.g., axisymmetry).
* **Other matter models** — we will include perfect fluids and vector fields in future.
---

## Summary of Features

| Feature                                  | Status          |
|------------------------------------------|-----------------|
| CCZ4 evolution                           | ✅ Ported       |
| Adaptive Mesh Refinement (AMR)           | ✅ Ported       |
| Checkpointing & Restart                  | ✅ Ported       |
| Moving puncture gauge                    | ✅ Ported       |
| BinaryBH example                         | ✅ Ported       |
| Diagnostic variables                     | ✅ Ported       |
| AMReX GPU offload                        | ✅ Ported       |
| MPI parallelism                          | ✅ Ported       |
| Scalar field matter model                | ✅ Ported       |
| Puncture tracking                        | ✅ Ported       |
| Particle-based diagnostics               | 🔧 In progress  |
| Weyl scalar / CCE extraction             | 🔧 In progress  |
| Apparent horizon finder                  | 🔧 In progress  |
| Initial condition solver                 | 🔧 In progress  |
| ADM quantities                           | 🔧 In progress  |
| 6th-order finite differencing            | 🔮 Planned      |
| Proca field matter                       | 🔮 Planned      |
| Perfect fluid matter                     | 🔮 Planned      |
| Cartoon formalism                        | 🔮 Planned      |

---

### Notes for Users

- GRTeclyn is designed as a research tool, not a black box. Users are encouraged to modify and extend it for their physics cases.
- If you're a student or supervisor considering using GRTeclyn, ensure your team has familiarity with numerical relativity, C++ and HPC environments.
- See the **Getting started** pages for how to compile and run your first simulation.

## License

GRTeclyn is made publicly available under the 3 clause BSD license.

### Copyright 2021, GRTL Collaboration

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Citations 
We would be grateful if, when using GRTeclyn, you cited our work. This ensures that all our developers get credit for the code. 

A paper for GRTeclyn will be coming soon, but in the meantime please include the citation for our GRChombo JOSS paper, which has the following bibtex reference:

```
@article{Andrade2021,
  doi = {10.21105/joss.03703},
  url = {https://doi.org/10.21105/joss.03703},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {68},
  pages = {3703},
  author = {Tomas Andrade and Llibert Areste Salo and Josu C. Aurrekoetxea and Jamie Bamber and Katy Clough and Robin Croft and Eloy de Jong and Amelia Drew and Alejandro Duran and Pedro G. Ferreira and Pau Figueras and Hal Finkel and Tiago Fran\c{c}a and Bo-Xuan Ge and Chenxia Gu and Thomas Helfer and Juha Jäykkä and Cristian Joana and Markus Kunesch and Kacper Kornet and Eugene A. Lim and Francesco Muia and Zainab Nazari and Miren Radia and Justin Ripley and Paul Shellard and Ulrich Sperhake and Dina Traykova and Saran Tunyasuvunakool and Zipeng Wang and James Y. Widdicombe and Kaze Wong},
  title = {GRChombo: An adaptable numerical relativity code for fundamental physics},
  journal = {Journal of Open Source Software}
}
```