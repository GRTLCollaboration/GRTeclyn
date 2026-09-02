# Building on GPUs

Compiling GRTeclyn for GPUs is in principle very easy, but in practise can be a bit of a pain. You will probably find it useful to look at some of the specific systems examples, e.g. the one for the Cosma8 AMD MI300 GPUs in [Example configs for specific HPC systems](example_configs.md).

The same process is followed as for CPUs, with the following changes:

## You need to install the right GPU compiler

One of the annoying things about GPUs is that there are 3 types related to the 3 vendors (Intel, AMD and Nvidia), and since they can't agree on things (thank you capitalism :pray:) you have to be able to understand the slightly different (but essentially the same) terminology and tools from all three.

Fortunately (thank you AMReX and historic US government funding :pray:), AMReX takes care of all the pain of implementing the code so it works on all three architectures. But you still have to think about this when you compile and run on a particular one. You may need to ask a lot of questions to the system admin, and you will probably end up feeling confused and stupid. That's ok. It will be worth it when you see the speed up, and once you are set up things should run smoothly.

For:

* AMD GPUs: you need to use the HIP backend
* Intel GPUs: you need to use SYCL backend
* Nvidia GPUs: you need to use the CUDA backend

You can think of them all as being like CUDA, but with a different name.

If you are on a login node, this will be a CPU node, but you may be able to module load the GPU compiler, e.g. on cosma8 you do
```
module load hipcc/6.3amd
```
If you are ssh-ing directly into the GPU node the compiler may be installed already by default, or you module load it as above. (TIP: You can always try to compile and see if it complains about not having the compiler.). For some computing clusters, it may be recommended to compile the code on the same GPU infrastructure that you intend on running on to prevent downstream issues. You may need to check with your system admins if this is required or not.

If you want to run with MPI over multiple GPUs, you will also need an MPI distribution, but note that this **needs to have been built with support for the GPUs you are using**, and often that won't be the case. You may need to ask the system admins for guidance. On Cosma8 they have an openmpi module with specific support, so you can just do:
```
module load openmpi/5.0.3
```

The good news about MPI is that if you really get stuck, it may be that you can just run on a single GPU, and don't need to use MPI. Incredibly, the BBH example runs on a single GPU really well. The problem with doing this is that usually the GPU node has multiple GPUs, and they are rarely set up for each one to support exclusive use. So you may find yourself sharing the memory of the GPU you are using with other users, which will utterly tank your performance. Therefore although in principle running on a single GPU without MPI is fine, in practise you may encounter problems until systems are configured better.

## You need to tell AMReX to compile with GPU offload support

This is the easy bit! Update your `Make.local-pre` to activate the appropriate options for the infrastructure you are using. For:

* AMD GPUs: set `USE_HIP=TRUE`
* Intel GPUs: set `USE_SYCL=TRUE`
* Nvidia GPUs: set `USE_CUDA=TRUE`

These options usually set a default compiler (`COMP`) which should work most of the time and there should rarely be a need to override this. Our [Example configs for specific HPC systems](example_configs.md) are a good place to look for inspiration.

You will also probably need to give it a specific flag for the GPU architecture, which you can google for (or hopefully find in the system docs). For example, your `Make.local-pre` may look something like:
```
USE_HIP=TRUE
# for AMD MI300
AMREX_AMD_ARCH=gfx942
# Optionally uncomment to turn off MPI
# USE_MPI=FALSE
```

Once you have created this file, you can then build in the usual way:
```
make -j 4
```

## You need to run it differently in the jobscript

Slurm isn't really designed for GPUs so the way you select the options in your jobscript can be a bit strange. Again it is worth asking for advice from the system admins if the documentation doesn't cover it, or looking at our example jobscripts. A usual set up is that you want to ask for one node that controls a certain number of GPUs, usually something like 8. As mentioned above, you would ideally like exclusive use of these GPUs, but if you pick a smaller number than the total number the node has that won't always be guaranteed.

The really important thing to understand is that how you are doing your parallelisation is very different now. GPUs are **huge and hungry** and they need to be fed a lot of points to process at the same time. So your grid is going to be divided up into a much smaller number of boxes, each with a lot of cells, and each big chunk will typically be given to a **single MPI process** running on a **single GPU**. This is why you can even use a single GPU to process the whole box in one go without using MPI at all. Make sure you have read [Performance optimisation](performance_optimisation.md) to understand how subdivision of the grid works, and consider whether you need to amend your params file to account for using GPUs (usually by increasing the max box size and blocking factor).

A typical command for running 8 MPI processes (which would be appropriate on a node that has 8 GPUs, and where you had 8 boxes to share out), is

```
mpirun -n 8 --exclusive ./main3d.hip.MPI.HIP.ex params.txt
```

Depending on the system you may also need to amend the `params.txt` file to specify the memory available to you and again make AMReX aware of the use of GPUs for how it uses MPI. e.g. for the MI300X I add:

```
# 192GB is the size of 1 GPU MI300X (it seems to work better to ask for a slightly smaller amount)
amrex.the_arena_init_size = 190000000000
amrex.use_gpu_aware_mpi = 1
```
(Usually you are using one MPI process per GPU, so you specify the per GPU amount of memory, not the total for your whole job.)
