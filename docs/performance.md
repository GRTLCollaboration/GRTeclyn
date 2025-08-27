## Load balancing on CPUs

GRTeclyn’s efficiency when running on a large number of distributed-memory nodes is highly dependent on good load balancing of the available computational work across those nodes. Load balancing seeks to avoid the situation where most of the nodes are waiting for some small subset of nodes to finish their computational work, and it does this by seeking to distribute the amount of work to be done per time step evenly among all of the nodes.

⚠️ Some of the load balancing work is done by AMReX automatically. 
**However, adjusting the input parameters will change how the main grids are split into boxes to be shared between processors, and understanding this and adapting the setup is crucial for the load balancing process.**

For example, if the coarsest level is divided into grid cells of 64^3, and you set the `max_box_size` to 64, then GRTeclyn *will not subdivide the coarsest grid*, but will just allow one box to cover the entire area. Thus only one process can work on this level, even if you run the code with 64 cores. If instead you set the maximum box size to 16, then the grid will be divided into (64/16)^3 = 64 boxes. Then running it on 64 cores, every process should get one box, and the problem will efficiently use the resources.

Consider also that running on 63 cores will mean that one process will act on two boxes, possibly taking twice as long to complete as the others. Since the processes must synchronise after each step, most will be waiting idly for this one to finish. So there is in principle no gain from using 63 cores rather than just 32. 

(Note that you might run on more cores because you need the additional memory, but you should still adjust the number of processes to match the load.)

Of course this calculation is much more complicated on the more refined levels where the number of boxes cannot necessarily be predicted ahead of time, but a bit of trial and error can still result in a big improvement. Note that the number of boxes on each process at leach level is output in the pout.x files, and so it is relatively easy to see how well load balanced you are by just running a few steps.

Note that load balancing the finest levels is *much* more important that balancing the lowest ones, since each finer level runs twice as often as the next coarser one.

There is also a minimum box size set by the `block_factor`, which we usually set to be equal (at least roughly) to the max box size, since this means that all the boxes are roughly the same size. Then having one box per process should mean roughly equal amounts of work. Below about `block_factor=8`, the costs of subdividing the grid start to outweigh the benefits of sharing the work (each box has +3 ghost cells on each edge, so the ghost cell load becomes comparable to the main calculation load). 

## Load balancing on GPUs

The really important thing to understand is that how you are doing your parallelisation is very different now. GPUs are **HUGE AND HUNGRY** and they need to be fed a lot of points to process at the same time. So your grid is going to be divided up into a much smaller number of boxes, each with a lot of cells, and each one will typically be given to an MPI process running on a single GPU. This is why you can even use a single GPU to process the whole box in one go without using MPI at all.
