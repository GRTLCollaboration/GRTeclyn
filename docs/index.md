# Welcome to GRTeclyn's documentation!

GRTeclyn is a flexible numerical relativity code that is the successor to [GRChombo](https://github.com/GRTLCollaboration/GRChombo), built on the [AMReX](https://amrex-codes.github.io/amrex/) library. The GRTeclyn source code is available on [GitHub](https://github.com/GRTLCollaboration/GRTeclyn) under the BSD-3-Clause license.

If you can't find the answers to your questions here, please post your questions on the [Discussions](https://github.com/GRTLCollaboration/GRTeclyn/discussions) section of the GitHub repo.

---

## Where should I start?

Here are some suggestions depending on your goal:

### *I want to know if GRTeclyn is the right code for my research project*

Start with the [**Capabilities**](capabilities.md) page for an overview of GRTeclyn's current features.

GRTeclyn is **not** a black-box code. You’ll need to understand how it works, and may need to write your own tools for specific applications. If you're a supervisor planning to assign a student to work on the code, please ensure you have sufficient background in **numerical relativity** and **C++/GPU computing** to support them.

For simpler projects e.g. for masters students, consider using the 1D code [**engrenage**](https://github.com/GRTLCollaboration/engrenage/wiki), which is a python based code that runs easily on laptops.


### *I want to compile and run an example (e.g. binary black hole merger)*

Check you have the [**Prerequisites**](prerequisites.md) and then follow the instructions for [**Building and running on CPUs**](building_cpus.md) or [**GPUs**](building_gpus.md). Our guides assume basic command-line experience (that of most new PhD students). We also assume access to a high-performance computing (HPC) cluster - while we provide laptop tips, most realistic simulations require HPC. 

Note that the code runs **a lot** faster on GPUs than on CPUs (up to a factor of 20), so if you can get access to a system with GPU nodes, this will be more efficient for most problems.

### *I already have NR experience and want to use or adapt the code*

1. Still begin with the [**Prerequisites**](prerequisites.md) and the [**Building and Running on CPUs**](building_cpus.md) or [**GPUs**](building_gpus.md) sections to ensure your environment is set up properly.

2. Then head to the section on [**Code Structure**](code_structure.md), which covers the design philosophy and guidance for adding your own physics.

For ex-GRChombo users, most of the code should look pretty familiar, but there is a helpful guide on [**Moving from GRChombo**](moving_from_grchombo.md) that you should start with.

If you have questions, feel free to [**Contact us**](contact_us.md) — we’ll try to help. But do keep in mind our unofficial motto:

> _“If you want it, you build it.”_

Before publishing results, don’t forget to check the [**License and citation**](license_and_citation.md) page.


### *I found a bug*

Please file an [issue](https://github.com/GRTLCollaboration/GRTeclyn/issues) on **GitHub**:

- Include your code version/commit hash
- Try to include as much information as you can, including: 
  - Detailed steps to reproduce the issue
  - Your build/runtime environment
  - Any error messages
  - The expected outcome/output of your run

GitHub Issues are preferred for code-related problems, but you can also post general questions on the [Discussions](https://github.com/GRTLCollaboration/GRTeclyn/discussions) section of the GitHub repo.


### *I want to contribute to GRTeclyn*

That's great! For now GRTeclyn is in development so it is not a good time for others to contribute, but we will post more details here once we are ready for your help.

---

## Disclaimer

**GRTeclyn is a research code under active development.**  
We aim to provide a stable core, but not all features are guaranteed to work perfectly.  
**Use it at your own risk!**
