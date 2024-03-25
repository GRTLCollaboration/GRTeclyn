This solves the Klein-Gordon equation $\frac{\partial^2 \phi}{\partial t^2} = \nabla \phi + V(\phi)$, based on the wave example in the [AMReX guided tutorials](https://github.com/AMReX-Codes/amrex-tutorials). 

The potential is $m^2 \phi^2$, and in the case of multiple scalar fields this is: $\phi^2 = \phi_1^2 + \phi_2^2 + ...$

The initial condition is a Gaussian pulse: $1.0 + A * exp(-r^2/\sigma)$


There are several properties of the scalar fields that are set in the input parameter file:
* wave.nfields - the number of scalar fields
* wave.initial_amplitude - the initial amplitude of the Gaussian pulse, there should nfield integers here, separated by a space
* wave.initial_width - intial width of Gaussian pulse, again one value per scalar field
* wave.scalar_mass - scalar field mass, one per scalar field





From the original README: 

>The Laplacian operator is discretized dimension by dimension with a fourth-order stencil,

>$$\frac{\partial^2 u}{\partial x^2} = \left(-\frac{5}{2} u_i + \frac{4}{3} (u_{i-1} + u_{i+1}) - \frac{1}{12} (u_{i-2} + u_{i+2})\right) / \Delta x^2$$

>The time stepping is done with a Runge-Kutta method (RK2, RK3 or RK4).  
>In this test, the displacement at the x-direction boundaries is zero, and the it's periodic in the y-direction.  
>Note that refluxing is not implemented in this test code.
