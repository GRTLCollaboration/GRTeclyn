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

* [Style guide for C++](https://google.github.io/styleguide/cppguide.html)
* [Basic guide to C++](https://www.tutorialspoint.com/cplusplus/index.htm)

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

<https://www.tutorialspoint.com/git/git_quick_guide.htm>

in particular, a useful guide to "rebase a branch to master" is

<https://www.theserverside.com/blog/Coffee-Talk-Java-News-Stories-and-Opinions/How-to-Git-rebase-a-branch-to-master-example>

## VisIt

The VisIt wiki is here:

<http://www.visitusers.org/index.php?title=Main_Page>

and the user manual is here:

<http://visit-sphinx-user-manual.readthedocs.io/en/latest/>

## yt

The yt package is also useful for analysing data, and has an excellent cookbook containing examples:

<http://yt-project.org/doc/index.html>

## Miscellaneous

A very useful GUI tool to compare two directories is Meld:

<https://meldmerge.org/>
