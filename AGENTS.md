# GRTeclyn Agents Guide

## Overview

GRTeclyn is a numerical relativity code for solving Einstein's equations using finite volume and adaptive mesh refinement methods. It is built on top of the [AMReX](https://github.com/AMReX-Codes/amrex) library but we use AMRLevel functions not the AMR functions directly. AMReX supplies the supporting functions for AMR, load balancing, the data containers, parallelization and GPU offloading.

## Purpose & Personas

- **GRTeclyn developers** - You are an expert in GPU programming. You follow a predictable, documented workflow that respects the coding conventions defined in the documentation in `docs/`. 
- **GRTeclyn users** - You will provide accurate build/test/docs guidance sourced from the repo.
- **Downstream module developers** — You are responsible for downstream repositories e.g. `GRFolres`, `GRDzhadzha`, that are built *on top of* GRTeclyn and pinned directly against its class/method names (`GRAmrLevel`, `SimulationParametersBase`, `StateVariables` enums, the `matter_t` interface `CCZ4RHSWithMatter` expects). As with GRTeclyn developers, but also be mindful of breaking GRTeclyn changes. 


## Documentation

Refer to the [website](https://grtlcollaboration.github.io/GRTeclyn/) for the code documentation. The website should be considered to be an authoritative source. 

## Code structure

The code is structured as follows:

- Source: This is where the main code is located. It is not specific to a particular physical scenario. 
- Tests: These are the unit tests run via [doctest](https://github.com/doctest/doctest). Place unit tests for new features here. 
- Examples: These are specific physics examples using the underlying code. BinaryBH runs a binary black hole merger. ScalarField runs an oscillaton. KleinGordon solves the Klein-Gordon equation without GR.
- Tools: Contains old doxygen documentation, the default top-level Makefiles called by the examples and the copyright info. 
- external: Contains the header for doctest. 
- docs: Documentation lives here. 
- .github: GitHub actions responsible for CI/CD. Also contains comparsion data for CI/CD.  
- .gitlab: GitLab pipeline that runs the regression tests on the A100 partition on CSD3. 


## Recommended workflow

When setting up new users to the code or on a new system/HPC cluster: 

- Install [pre-commit](https://pre-commit.com/) including clang-format and clang-tidy. Refer to `.pre-commit-config.yaml` for the correct version number of both.
- Git clone AMReX 
- Git clone GRTeclyn. Prefer the GRTeclyn directory to live next to the AMReX directory not as a sub-directory or vice versa.
- Make a new branch from `main` for changes.
- Refer to the coding conventions in the documentation.
- Manually run the unit and regression tests after each commit. 

## Rules for pull requests

- Never push directly to `main`. This is protected by branch protection rules also. 
- New feature branches should be named `enhancement/my-feature`.
- Bugfixes should be named `bugfix/github-issue-number-or-name`. If no relevant issue number exists, create an issue.
- Use a descriptive title with less than 50 characters.
- Write a descriptive message that is formatted to have less than 70 characters per line.
- Require that all the CI/CD tests have passed, except for draft pull requests. 
- Before opening a pull request, git pull from the remote to ensure that the branch is up to date. Rebase if necessary. 
- Small commits e.g. changes for the linter, should be squashed or fixed up.
- Remind the user to update the documentation with the pull request. Always update the documentation if changes made by the pull request also update information in the documentation. Always update the documentation if changes made by the pull request are breaking changes.


## Rules for opening a new issue

- Check if the branch is up to date with `main`.
- Check if there is an existing open issue that covers the same topic.
- Write a detailed report on how to reproduce the issue, including the build environment and flags used, and quote the GitHub SHA.


## GPU Lambda Safety

Identical caveat to AMReX's own `AGENTS.md`

