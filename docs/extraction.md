# Extraction

## Overview

Before reading this section on extracting various fields using a specified geometry, we recommend you check out the [Particle Interpolator](https://grtlcollaboration.github.io/GRTeclyn/particle_interpolator/) chapter. If you would like to use gravitational-wave (GW) extraction or perhaps adapt the current extraction classes to your own needs, this chapter is for you.

!!! note

    Extraction uses our `ParticleInterpolator` class. Hence, if we want persistent extraction throughout the simulation, we should really initiate the `ParticleInterpolator` only once. An example of this is explained in GW extraction section below.

We inherit our legacy GRChombo structure, where we have a general class called `SurfaceExtraction`, which sets up the prerequisites for a 2D surface extraction. The geometry and the variables you are extracting in this class are kept arbitrary, but we will explain how they need to be specified shortly. Function `SurfaceExtraction::extract()` contains the necessary call to the `ParticleInterpolator` class; recall that your extraction points may not coincide with the grid points. We note that the class takes the following form:

```cpp
template <class SurfaceGeometry, int num_components> class SurfaceExtraction
```

So, similar to GRChombo, we template over the `SurfaceGeometry`, which can be e.g. `SphericalGeometry`, but also add an additional template `num_components`, which directly enters the `ParticleInterpolator` and corresponds to the number of components we are extracting. The `SurfaceExtraction` class additionally has several helper functions that let you add the variables you need to extract, see e.g. `SurfaceExtraction::add_state_vars()` and `SurfaceExtraction::add_derived_vars()`. The information about the variables that you are extracting is stored in a special struct:

```cpp
struct vars_t
    {
        int comp{};
        VariableType type{};
        Derivative deriv;

        amrex::Vector<BCParity> parities{};
        std::string derived_name;
    };
```

which then propagates all of this information directly to the `ParticleInterpolator`. Why is this needed? Well, you should do your homework and read the [Particle Interpolator](https://grtlcollaboration.github.io/GRTeclyn/particle_interpolator/) chapter :zany_face: !

We similarly reuse GRChombo's `SphericalExtraction` class, which allows for the spin-weighted spherical-harmonic decomposition. Similar to `SurfaceExtraction` it also now has an additional template argument for `int num_components`.

## GW extraction

The most relevant application of the extraction classes to numerical relativity is undoubtedly the GW extraction. Let us walk you through how we implement this, as it may be useful if you wish to extract other quantities.

### Specifying how to extract $\Psi_4$

The relevant class to start with is `WeylExtraction`, where provide the machinery to extract $\Psi_4$.

First, our class is defined as follows:

```cpp
class WeylExtraction : public SphericalExtraction<2>
```

Here:

* we inherit from `SphericalExtraction`, since we extract on spheres. This fixes the geometry to `SphericalGeometry` and gives us access to additional post-processing functions for mode decomposition.
* `SphericalExtraction` is templated over two components. Recall that $\Psi_4$ has two components -- real and imaginary.

Within the class, we add components needed to be extracted and note that here we additionally provide the parities. This is important because $\Psi_4$ is a **derived** variable, not a state variable!

```cpp
amrex::Vector<BCParity> parities = {BCParity::even, BCParity::odd_xyz};
this->add_derived_vars({0, 1}, parities, "Weyl4");
```

Here `{0, 1}` refers to the component indices, whilst `Weyl4` refers to the derived group. We encourage you to consult the `Weyl4` class for this.

!!! warning

    Our extraction implementation is such that you cannot mix different derived groups. You will need to define separate objects if you wish to extract several derived groups!


We then provide the instruction to the class on what to do, i.e. how we want to extract $\Psi_4$ using `WeylExtraction::execute_query()`. In particular,

* we extract $\Psi_4$ on the sphere using `extract(a_interpolator)`. This can also write out the field values at the specified $(r, \theta, \phi)$ coordinates for each time step if you set `write_extraction` to `true`.
* we additionally perform mode decomposition using `add_mode_integrand()` and `integrate()`.

### Linking with your example and `ParticleInterpolator`

In this section we will explain how to link your extraction method to your particlar example. We will use the BinaryBH as an example.

As already noted in the blue box at the beginning of this chapter, we should ideally instantiate the `ParticleInterpolator` only once. We therefore store it as part of the AMR object and set it up during AMR initialisation. For example, for black holes (BHs) we have `BHAmr` class. Hence, we directly add the interpolator as its public member:

```cpp
// example for interpolator object for Psi4 extraction
    static constexpr int weyl_num_components = 2;
    ParticleInterpolator<weyl_num_components> m_weyl_interpolator; // interpolator object
```

and set up the interpolator in `BHAmr::init()`:

```cpp
void init(amrex::Real a_strt_time, amrex::Real a_stop_time) override
    {
        GRAmr::init(a_strt_time, a_stop_time);

        const auto &params = get_simulation_parameters();
        m_weyl_interpolator.setup(this, params.boundary_params, true);
    }
```

In our binary BH example, we can then switch on the extraction in `BinaryBHLevel::specific_post_timestep` using:

```cpp
// Weyl extraction
    if (extraction_params.enabled)
    {
        int min_level = extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);

        if (calculate_weyl && Level() == min_level)
        {
            amrex::Real m_time       = get_state_data(state_index).curTime();
            amrex::Real m_dt         = get_gr_amr_ptr()->dtLevel(Level());
            amrex::Real restart_time = get_gr_amr_ptr()->get_restart_time();
            bool first_step          = (m_time <= m_dt);

            WeylExtraction my_extraction(extraction_params, m_dt,
                                         m_time, first_step, restart_time);
            my_extraction.execute_query(&get_bhamr_ptr()->m_weyl_interpolator);
        }
    }
```

!!! note

    If you are using derived variables, such as $\Psi_4$, do not forget to register the derived quantity in your `Example::variableSetUp()`. For instance, for $\Psi_4$, one would add `Weyl4::set_up(state_index)`.

## Defining your own extraction

It is hard to cover every single extraction scenario a user may want to have, as it depends on the physics one is trying to do. Here we briefly describe what you need to do in the most common situations:

1. *I want to extract my own quantity on the sphere.*

    Please refer to the `WeylExtraction` class and the previous section for an example that you can adapt to your needs.

2. *I want to use my own geometry/extraction.*

    If your application requires something other than spheres, you should write your own `MyGeometry` class and `MyExtraction` class. In fact, our GRChombo code has an example of [cylindrical extraction](https://github.com/GRTLCollaboration/GRChombo/blob/training/fawcett/Source/AMRInterpolator/CylindricalExtraction.hpp) and the associated [cylindrical geometry](https://github.com/GRTLCollaboration/GRChombo/blob/training/fawcett/Source/AMRInterpolator/CylindricalGeometry.hpp). However, this has not been ported to GRTeclyn yet.

3. *Can I extract along the line?*

    Yes! We provide a `LineExtraction` class that does exactly that and requires the user to specify the start and end points of the line. This is also a nice example to see how `ParticleInterpolator` is used.





