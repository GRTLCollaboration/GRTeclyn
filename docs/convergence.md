> “Whatever the ways of the universe may require us to suffer, let us take it up with high mindedness. This is the oath by which we are bound, not to be disturbed by what we do not have the power to avoid.” Seneca, On the Happy Life

Seneca was certainly thinking about convergence testing when he said this.

Convergence testing with an AMR code can be challenging, but it is essential. The basic principle of convergence testing is that, once the solution is well resolved, running with a higher resolution should give roughly the same result, with an error that decreases at an order set by the finite difference and time evolution stencils used, in the case of GRTeclyn this is 4th order.

Here is a rough guide to how we do convergence testing. You should refer to Alcubierre's book on numerical relativity for a more detailed treatment of convergence and convergence rates. Here we just give practical tips related to GRTeclyn.

- *What should I do my convergence testing on?*

  The general rule is that you should convergence test all the physical results that you want to present in any research publication for which you are using the code. So if you are calculating GW waveforms, you need to check their amplitude and phase converge. If you are tracking a total mass of matter, you should check that. In general global quantities will give better convergence than local ones, but if you are showing plots of the field value at a particular point, you might want to convergence test that.

  One thing you can and should always do a convergence test on is the constraint error - this has the advantage that you know it should be zero, so you can calculate the convergence rate with only 2 (rather than 3) simulation resolutions. Usually we take a norm of the constraint across the region of interest - this may involve zeroing it near the boundaries, where you will always get constraint error because Sommerfeld conditions are not constraint preserving. This may feel like cheating, but as long as the region of physical interest is converging well, this is good enough to support your physical result. Likewise it is conventional to zero within any BH horizons to remove the spurious effect of the singularity (which has higher constraint error the better your resolution).

  Note that there are functions for taking quantities across the grid in the AMRReductions class. These will often be very useful for convergence testing as they allow global quantities to be extracted.

- *How do I set the resolutions so the grid does the right thing?*

  This is where AMR makes thing hard. Say you are trying to double the resolution (normally in practise you would aim for some smaller ratio like 1.5 as doubling is expensive, but let's say double just for illustration). Your goal should be to double the resolution of every region of the domain, so each grid should cover the same area as before but with double the resolution. Obviously you need to double the resolution of the coarsest grid, and reduce your tagging thresholds, but depending on the tagging criteria you are using getting regions to match up may be very hard to achieve exactly. However, provided you roughly cover the same regions at double resolution, and in particular make sure any regions driving your convergence quantity are appropriately re-refined, the convergence test should work ok. In general we aim to adjust the thresholds until we get the finest grid to cover the same area at double the resolution, and let the other levels do what they want to do.

- *I have run three versions of the simulation at 3 different resolutions HR, MR and LR, now what?*

  You now need to take your outputs at each resolution and difference them, so ErrorH = HR - MR and ErrorL = MR - LR. ErrorL should be bigger than ErrorH - if it is not you are in trouble and definitely need higher resolution. Note that you may need to use an interpolating function (e.g. in Python) to obtain comparable data at the same times, since the time steps won't match up at different resolutions.

  However, just seeing that error decreases with resolution *tells you nothing*. The expected order of the convergence must be right - so for 4th order stencils and doubling the resolutions it should be (1/2)^4 smaller.   This validates that this specific source of error (ie, from the numerical integration) is dominating, which is the goal of the convergence test.   If you obtain a different rate, some other error source is dominating, and if you don’t know what that is, you could be in trouble.   

  e.g. If you hit a place where you converge slower than expected, it could be that numerical precision error is now dominating. This might be ok, but then it also tells you that you are wasting time resolution, because you have hit a wall where you can’t get better by taking better resolutions.

  Equally, if you find the error getting smaller too fast, don't celebrate - it is likely that your low resolution run was not in the convergence regime, so you should again increase resolution to check you can achieve the appropriate rate.

- *I don't get exactly the right convergence rate, but my results look exactly the same at each resolution and when I add more resolution I get the same thing. Help!*

  This is when you need to apply some judgement. Anyone who does simulations will tell you that convergence is an art rather than a science, and it won't always be perfect. If you plot your quantity of interest and you can't see any difference in the 3 resolutions by eye (and you have tried increasing the resolution again with the same result), it might be fine and you have just reached some fundamental limit of numerical precision somewhere. But higher resolution should always give you smaller error, otherwise you definitely need to investigate more. You might especially lose some convergence order when there is a lot of rapid regridding in the AMR, which does introduce errors.

  You should also consider the accuracy of your initial data - if you use data from the [GRTresna Initial Condition Solver](https://github.com/GRTLCollaboration/GRTresna) that is only 2nd order accurate, and this error will tend to dominate at initial times, so you may only recover 3rd/4th order at later times once the evolution error starts to dominate.

  Also consider whether you have used some interpolation error in extracting the data or taking a volume integral - this will add a source of error too. The AMRInterpolator uses by default 4th order interpolation, so should not spoil the convergence results. However, yt uses zeroth order interpolation unless you tell it to do otherwise.