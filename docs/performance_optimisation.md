# Performance optimisation

GRTeclyn’s efficiency when running on a large number of distributed-memory nodes is highly dependent on good load balancing of the available computational work across those nodes. Load balancing seeks to avoid the situation where most of the nodes are waiting for some small subset of nodes to finish their computational work, and it does this by seeking to distribute the amount of work to be done per time step evenly among all of the nodes.

Some of the load balancing work is done by AMReX automatically.
**However, adjusting the input parameters will change how the main grids are split into boxes to be shared between processors, and understanding this and adapting the setup is crucial for the load balancing process.**

## Load balancing on CPUs

For example, if the coarsest level is divided into grid cells of 64^3, and you set the `max_box_size` to 64, then GRTeclyn *will not subdivide the coarsest grid*, but will just allow one box to cover the entire area. Thus only one process can work on this level, even if you run the code with 64 cores. If instead you set the maximum box size to 16, then the grid will be divided into (64/16)^3 = 64 boxes. Then running it on 64 cores, every process should get one box, and the problem will efficiently use the resources.

Consider also that running on 63 cores will mean that one process will act on two boxes, possibly taking twice as long to complete as the others. Since the processes must synchronise after each step, most will be waiting idly for this one to finish. So there is in principle no gain from using 63 cores rather than just 32.

(Note that you might run on more cores because you need the additional memory, but you should still adjust the number of processes to match the load.)

Of course this calculation is much more complicated on the more refined levels where the number of boxes cannot necessarily be predicted ahead of time, but a bit of trial and error can still result in a big improvement. Note that the number of boxes on each process at leach level is output in the pout.x files, and so it is relatively easy to see how well load balanced you are by just running a few steps.

Note that load balancing the finest levels is *much* more important that balancing the lowest ones, since each finer level runs twice as often as the next coarser one.

There is also a minimum box size set by the `block_factor`, which we usually set to be equal (at least roughly) to the max box size, since this means that all the boxes are roughly the same size. Then having one box per process should mean roughly equal amounts of work. Below about `block_factor=8`, the costs of subdividing the grid start to outweigh the benefits of sharing the work (each box has +3 ghost cells on each edge, so the ghost cell load becomes comparable to the main calculation load).

## Load balancing on GPUs

The really important thing to understand is that how you are doing your parallelisation is very different now. GPUs are **HUGE AND HUNGRY** and they need to be fed a lot of points to process at the same time. So your grid is going to be divided up into a much smaller number of boxes, each with a lot of cells, and each one will typically be given to an MPI process running on a single GPU. This is why you can even use a single GPU to process the whole box in one go without using MPI at all.

## Writing code for GPUs

There are some factors that should be considered when writing GPU performant code. 

### Limit the amount of memory used

This seems like a simple one, but having a large memory footprint affects GPUs differently to CPUs. The GPUs contain contain simpler processors and do not have an ability to "pre-fetch" instructions unlike CPUs. 

Just like CPUs, the GPUs have a memory hierarchy. However, unlike (most) CPUs, local memory is not a software abstraction. Local memory to each thread is a dedicated resource and is the fastest to read and write, follwed by shared work group/thread block memory then global memory (also called VRAM). But local memory also has the smallest volume. If you reduce the amount of memory needed in your calculation, you can encourage the compiler to keep your variables stored on the *lowest* level of memory possible. This is why we no longer have the `Vars` struct in GRChombo or the derivative structs - these were being stored and fetched inside global memory and slowing down the code.

!!! warning

    Only instantiate the variables that you are using at the time that they are needed. Limit the persistence of variables.


### Kernel Fission

Related to the issue above, if the calculation uses too much local (register) memory, memory spills can occur. You can see this at compile time, e.g. for CUDA, there will be a line reporting the number of bytes spilled. 

This is why the `CCZ4RHS` calculation has been split up into three different parts:

- one `ParallelFor` launch to calculate `chi` and `h_ij`
- one `ParallelFor` launch to calculate `A_ij`, `Theta` and `Gamma_ij`
- one `ParallelFor` launch for the gauge and dissipation

Take a look inside `BinaryBHLevel.cpp` in the `specificEvalRHS` member function as an example of how to do the splitting. 

!!! warning

    Use kernel splitting when calculating the RHS. Always call the above functions in separate `ParallelFor` instances. 


### Kernel Fusion

In circumstances where there is not enough work for the GPU to do, your code may benefit from the opposite effect, kernel *fusion*. Small workloads, such as initializing variables, should be combined together in the same `ParallelFor` launch whenever possible. This reduces the launch overhead from each call to `ParallelFor`. 

!!! warning

    Try to combine small workloads into a single `ParallelFor` launch.


### Tensors 

The `Tensor` object has been overhauled since GRChombo. It is now a lightweight wrapper around `amrex::Array1D`, `amrex::Array2D` or `amrex::Array3D` depending on the rank of the `Tensor` object.

It supports curly brace initialization and you can access the elements using `()` brackets, e.g.:
```
Tensor::Rank2 my_tensor{0.};     // all values are initialized to 0.
FOR(i)
{
    my_tensor(i, i) = 1;         // set the diagonal elements to 1. 
}
```


The default size of the tensor object is `AMREX_SPACEDIM` but you can have any size that you like if you write it like this:

```
// The syntax is Tensor::GeneralRank<RANK, DIM1, DIM2, DIM3>

Tensor::GeneralRank<3, 10, 20, 30> my_oddly_shaped_tensor{0.}; 

my_oddly_shaped_tensor(1, 2, 3) = 4. // still OK!
```

You can declare tensors up to rank 4, after that you will need to modify `Tensor.hpp` to define your own.

There is even a special rank-4 tensor with size `AMREX_SPACEDIM+1`:

```
Tensor::SpacetimeRank4 my_spacetime_tensor{0.}

for (int i = 0; i < AMREX_SPACEDIM+1; ++i)   // NB: can't use FOR!
  for (int j = 0; j < AMREX_SPACEDIM+1; ++j)
    for (int k = 0; k < AMREX_SPACEDIM+1; ++k)
      for (int l = 0; l < AMREX_SPACEDIM+1; ++l)      
        my_spacetime_tensor(i, j, k, l) = i + j + k + l;  // or something...

```

Some variables are symmetric so we only store the unique values. This is the same as GRChombo where `Vars` only stored `A12` but not `A21`. For example: 

```
Tensor::Sym12Sym34Rank4 my_symmetric_tensor{0.}; //initialize all elements to 0

FOR(i)
 FOR(j)
  FOR(k)
   FOR(l)
    my_symmetric_tensor(i, j, k, l) = i + cos(j) + log(k) + l;  

```
In the example above, not all the values will be saved to `my_symmetric_tensor` as the calculation is not symmetric in the `ij` or `kl` pair of indices. 

Some things to bear in mind: 

- The second derivatives are always assumed to be symmetric: $\frac{d^2}{dxdy} = \frac{d^2}{dydx}$ and we will only store the unique values. 
- `h_ij`, `d1_h` and `d2_h` are always symmetric. `d2_h` is special in that it is doubly symmetric so it has type `Tensor::Sym12Sym34Rank4`. `d1_h` is only symmetric in $h_{ij}$ and not the derivative part, so `d1_h` has type `Tensor::Sym12Rank3` because the derivative is always the last index.
- `A_ij` and `d1_A` are always symmetric, similarly `d1_A` has type `Tensor::Sym12Rank3`. 
- The real dimension of `Tensor::Sym12Sym34Rank4` is 2. It is actually an 6 x 6 `amrex::Array2D` object (which is actually a 1D object...) as there only 36 unique values once both symmetries have been accounted. But you index into it the same as a 4D object, so always write `d2_h(i, j, k, l)` the code will work out the correct element. 


!!! warning

    Take advantage of symmetries whenever possible. 


### Limit branching within a kernel

  Because of their simpler processors, GPUs have limited branch prediction capabilities. Control flow statements, such as `if` statements, can severely impact performance by reducing thread occupation, as some branches are masked off or predicated in a different instruction stream. In CUDA terms, this is called "warp divergence".

If you must use branching, try to reword it as `if constexpr my_expression` and have `my_expression` be evaluated at compile time instead of runtime. 

!!! warning

    Try not to use `if...else...` statements inside a `ParallelFor`. Prefer the use of `if constexpr ...`


### Further reading

- **Data Parallel C++: Mastering DPC++ for Programming of Hetrogeneous Systems using C++ and SYCL**, Reinders, Ashbaugh, Brodman, Kinsner, Pennycook, Tian. A great intermediate level book, not limited to SYCL, one of Juliana's favourites. [Download](https://link.springer.com/book/10.1007/978-1-4842-5574-2)

- **Software Carpentries Lesson** Introductory material for GPU programming. If this is the first time you have ever programmed for GPUs, this is for you. [Start here](https://carpentries-incubator.github.io/lesson-gpu-programming/introduction.html)