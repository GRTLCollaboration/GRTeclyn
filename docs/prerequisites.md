To compile and run **GRTeclyn**, you will need the following software and tools:

---

## 🧰 Required Software

- **C++ compiler with C++17 support**
  GRTeclyn uses modern C++17 features. Supported compilers include:
  - GCC ≥ 8
  - Clang ≥ 6
  - Intel Classic ≥ 19.1.4

⚠️ The use of C++17 is stricter than GRChombo’s C++14 requirement. Older compilers that worked with GRChombo may not work here — but if you have a newer one available, and could build GRChombo, you can likely build GRTeclyn too.

- **[MPI](https://www.mpi-forum.org/)** (optional but recommended)
  Required for running in parallel across nodes or cores.
  Supported implementations: OpenMPI, MPICH, Intel MPI, etc.

Many computer clusters use some form of `modules` command to control the environment. They will probably have the C++ and MPI prerequisites installed above but you will need to manually load them (and unload any preloaded conflicting ones).

Here are some useful module commands:

* `module avail [<name>]`
  * What modules are available/available with name `<name>`?
* `module [un]load <module name>`
  * Load/unload a module.
* `module swap <old module> <new module>`
  * Swap `<old module>` for `<new module>`.
* `module list`
  * What modules are loaded currently?
* `module display <module name>`
  * Tell me more about `<module name>`: how does it change the environment?

Some newer clusters use [spack](https://spack.io/) to manage the installation of modules. On these systems, modules installed with spack will typically be appended by a 7 character hash (e.g. `hdf5-1.10.4-gcc-5.4.0-7zl2gou` on the CSD3 cluster). If there are multiple modules with similar names that only differ by their hash, you can can query the dependencies with the command `spack find -lvd <name>`.

> **Warning**  See [this page](https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#gcc-on-macos) if you want to build on macOS using Homebrew GCC.

The following are required, but are likely to be loaded by default on most systems:

- **[Git](https://git-scm.com/)**
  Used to obtain GRTeclyn and AMReX.

- **[GNU Make](https://www.gnu.org/software/make/)** (version ≥ 3.81)
  Standard build system used across the codebase.

- **[Python](https://www.python.org/)** (version ≥ 2.7)
  Used in build and analysis utilities.

- **Unix-like environment** with:
  - `perl`
  - `sed`
  - Standard shell tools
  These are available by default on most Linux/macOS systems.

---

## 🛠️ Obtaining the AMReX and GRTeclyn codes

First `cd` into a directory you are happy to clone the code into. For simplicity, we shall assume that is your home directory (so adjust any commands below accordingly if not).

The AMReX source code is hosted on [GitHub](https://github.com/AMReX-Codes/amrex). Clone this with the command

```bash
git clone https://github.com/AMReX-Codes/amrex.git
```
It will be cloned to the `amrex` directory.

Next, clone the GRTeclyn repository
```bash
git clone https://github.com/GRChombo/GRTeclyn.git
```

**Note**
We have assumed that you have cloned both of these repositories to the same directory so that the `amrex` and `GRTeclyn` directories share the same parent directory. This fact is used to locate AMReX when GRTeclyn builds its examples, and it is recommended as it makes it easier to be sure about the version you are using.

If you want to keep AMReX elsewhere, you need to set the `AMREX_HOME` environment variable appropriately in your shell
```bash
export AMREX_HOME=/path/to/amrex
```
(and you probably want to add this to your `.bashrc` file in your home directory too so it is run automatically when you log into the cluster you are using).

---

## Setup the pre-commit hook
We like to maintain a certain format to our coding style for uniformity and correctness which is enforced by `clang-format`. This can be checked within your editor and should be checked again when committing new code to the repository. Please see this [page](https://github.com/GRTLCollaboration/GRTeclyn/wiki/Using-Clang-Format) for more details on how to install `clang-format` and 'pre-commit'.

The pre-commit hook protects us from writing badly formatted code and should be installed whenever you are working on a new system.