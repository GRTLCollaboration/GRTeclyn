# Some examples 

## Binary Black Hole merger

This is the classic binary black hole merger from GRChombo.

Start here if you:

- want to learn how to port your code from GRChombo
- need to evolve all CCZ4 variables


## A simple Klein Gordon solver

This example solves the Klein Gordon equation without GR.

Start here if you:

- find the CCZ4 equations confusing and just want to learn the structure of the code and AMReX
- have your own RHS that you would like to implement separately for several classes
- want to work on optimizating the code (the register pressure is too high on the PVCs for an analysis using the BinaryBH example) 