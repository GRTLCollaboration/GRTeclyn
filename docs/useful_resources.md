# Useful resources

## AMReX

GRTeclyn is built on top of the AMReX library and we try to follow the thinking behind AMReX pretty closely e.g. the GPU strategy. 
The [AMReX documentation](https://amrex-codes.github.io/amrex/docs_html/) is a great resource that introduces a lot of basic concepts that you will encounter in GRTeclyn e.g. `MultiFab`, `Array4` and `ParallelFor` as well as build options for various GPU architectures. There are also suggestions for analysing and visualising PlotFiles and profiling and debugging. If a problem you encounter is not covered in this (admittedly short) wiki, then the AMReX documentation should be your next stop. 

There is also a set of short demos and exercises compiled by the AMReX team available [here](https://github.com/AMReX-Codes/amrex-tutorials).

Failing that, the tests directory of AMReX, `${AMREX_HOME}/Tests`, also functions as code examples for how to use certain AMReX features. 

## Numerical Relativity (NR)

There are many useful introductions to NR, but we recommend in particular the following open source resources:

* [Living Reviews in Relativity](http://www.springer.com/gp/livingreviews/relativity/lrr-articles) has a dedicated section on Numerical Relativity.
* [Eric Gourgoulhon's notes on NR](https://arxiv.org/abs/gr-qc/0703035) give a detailed introduction to the ADM decomposition, with a focus on the geometric interpretation.

The textbooks of Alcubierre, Baumgarte & Shapiro, and Shibata are also bookcase essentials for anyone working in NR.

## Useful C++ resources

Style guide for C++
https://google.github.io/styleguide/cppguide.html

Basic guide to C++
https://www.tutorialspoint.com/cplusplus/index.htm

One should review in particular the sections on classes and templating, which are both used extensively in the code.


All these places offer courses relevant for scientific programming, best check them regularly 

* https://www.archer.ac.uk/training/ 
* https://training.prace-ri.eu/
* https://extremecomputingtraining.anl.gov/ ( ATPESC - yearly course, typically around a week. Recommended for more advanced students)
* https://learn.tacc.utexas.edu/ 
* https://portal.xsede.org/training/overview 

## HPC cluster usage

Most machines have their own webpages which give a guide to usage, and tell you things particular to each machine, like compilers available, the job submission system, where to store your data etc.

The Rosetta stone of job submission commands is here:
https://slurm.schedmd.com/rosetta.pdf

Some useful commands:

    - 'passwd' can be used to change your login password.
    - 'grep' to search for text in your files 

## Git guide

A comprehensive guide to git:

https://www.tutorialspoint.com/git/git_quick_guide.htm

in particular, a useful guide to "rebase a branch to master" is

https://www.theserverside.com/blog/Coffee-Talk-Java-News-Stories-and-Opinions/How-to-Git-rebase-a-branch-to-master-example

## VisIt

The VisIt wiki is here:

http://www.visitusers.org/index.php?title=Main_Page

and the user manual is here:

http://visit-sphinx-user-manual.readthedocs.io/en/latest/

## yt

The yt package is also useful for analysing data, and has an excellent cookbook containing examples:

http://yt-project.org/doc/index.html

## Miscellaneous

A very useful GUI tool to compare two directories is Meld:

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
https://meldmerge.org/
