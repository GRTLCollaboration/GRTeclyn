/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SURFACEEXTRACTION_HPP_
#define SURFACEEXTRACTION_HPP_

// AMReX includes
#include <AMReX_GpuContainers.H>

// Parameters
#include "SurfaceExtractionParameters.hpp"

// Other includes
#include "DimensionDefinitions.hpp"
#include "FilesystemTools.hpp"
#include "IntegrationMethod.hpp"
#include "InterpolationQueryParticle.hpp"
#include "Lagrange.hpp"
#include "ParticleInterpolator.hpp"
#include "SmallDataIO.hpp" // for writing data
#include "StateTypes.hpp"
#include "StateVariables.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <utility>
#include <vector>

//! This class extracts grid variables on 2 dimensional surfaces each
//! parameterised by u and v with different surfaces given by level sets of
//! another parameter
template <class SurfaceGeometry, int num_components> class SurfaceExtraction
{
  public:
    // I suggest using a struct here instead of the tuple-structure inherited
    // from GRChombo here, as we need to add more features to the varibles used
    struct vars_t
    {
        int comp{};
        VariableType type{};
        Derivative deriv;

        amrex::Vector<BCParity> parities{};
        std::string derived_name;
    };
    using params_t = surface_extraction_params_t;

  protected:
    SurfaceGeometry m_geom; //!< the geometry class which knows about
                            //!< the particular surface
    params_t m_params;
    std::vector<vars_t>
        m_vars; //!< the vector of of variables and their features to extract
    amrex::Real m_dt{};
    amrex::Real m_time{};
    bool m_first_step{};
    amrex::Real m_restart_time{};
    int m_num_interp_points{};  //!< the total number of points this
                                //!< rank will extract (0 on ranks > 0)
    amrex::ParticleReal m_du{}; //!< the grid spacing in u (used in integrate)
    amrex::ParticleReal m_dv{}; //!< the grid spacing in v (used in integrate)

    std::vector<std::vector<amrex::ParticleReal>> m_interp_data;
    std::array<std::vector<amrex::ParticleReal>, AMREX_SPACEDIM>
        m_interp_coords;
    // this is the really long type used for integrands
    // the vector<amrex::ParticleReal> is a vector of all the extracted
    // variables at that point in the order they were added
    using integrand_t = std::function<amrex::ParticleReal(
        std::vector<amrex::ParticleReal> &, amrex::ParticleReal,
        amrex::ParticleReal, amrex::ParticleReal)>;
    std::vector<integrand_t> m_integrands;
    std::vector<std::array<IntegrationMethod, 2>> m_integration_methods;
    std::vector<std::reference_wrapper<std::vector<amrex::ParticleReal>>>
        m_integrals;
    std::vector<bool> m_broadcast_integrals;

    bool m_done_extraction{}; //!< whether or not the extract function has
                              //!< been called for this object

    //! returns the flattened index for m_interp_data and m_interp_coords
    //! associated to given surface, u and v indices
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    [[nodiscard]] int index(int a_isurface, int a_iu, int a_iv) const
    {
        return a_isurface * m_params.num_points_u * m_params.num_points_v +
               a_iu * m_params.num_points_v + a_iv;
    }

  public:
    //! Normal constructor which requires vars to be added after construction
    //! using add_var or add_vars
    SurfaceExtraction(surface_extraction_params_t a_params, amrex::Real a_dt,
                      amrex::Real a_time, bool a_first_step,
                      amrex::Real a_restart_time = 0.0);

    //! Constructor for geometries which require instance-specific parameters
    SurfaceExtraction(surface_extraction_params_t a_params,
                      SurfaceGeometry a_geom, amrex::Real a_dt,
                      amrex::Real a_time, bool a_first_step,
                      amrex::Real a_restart_time = 0.0);

    //! add a single variable or derivative of variable
    void add_var(int a_var, const VariableType var_type = VariableType::state,
                 const Derivative &a_deriv                 = Derivative::LOCAL,
                 const amrex::Vector<BCParity> &a_parities = {},
                 const std::string &a_derived_name         = "");

    //! add a vector of variables/derivatives of variables
    void add_vars(const std::vector<vars_t> &a_vars);

    //! add a vector of state variables (no derivatives)
    void add_state_vars(const std::vector<int> &a_vars);

    //! add a vector of derived variables (no derivatives)
    void add_derived_vars(const std::vector<int> &a_vars,
                          const amrex::Vector<BCParity> &a_parities,
                          const std::string &a_name_derived);

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    //! Alternative constructor with a predefined vector of variables and
    //! derivatives
    SurfaceExtraction(const params_t &a_params,
                      const std::vector<vars_t> &a_vars, amrex::Real a_dt,
                      amrex::Real a_time, bool a_first_step,
                      amrex::Real a_restart_time = 0.0);

    //! Another alternative constructor with a predefined vector of variables
    //! no derivatives
    SurfaceExtraction(const params_t &a_params, const std::vector<int> &a_vars,
                      amrex::Real a_dt, amrex::Real a_time, bool a_first_step,
                      amrex::Real a_restart_time = 0.0);
    // NOLINTEND(bugprone-easily-swappable-parameters)

    //! Do the extraction
    void extract(ParticleInterpolator<num_components> *a_interpolator);

    //! Add an integrand dependent on the interpolated data over the surface
    //! for integrate() to integrate over.
    //! Note the area_element is already included from the SurfaceGeometry
    //! template class
    //! The last argument is whether to broadcast the result to all MPI ranks
    //! or just keep on rank 0. Most use cases won't need this set to true.
    void add_integrand(
        const integrand_t &a_integrand,
        std::vector<amrex::ParticleReal> &out_integrals,
        const IntegrationMethod &a_method_u = IntegrationMethod::trapezium,
        const IntegrationMethod &a_method_v = IntegrationMethod::trapezium,
        const bool a_broadcast_integral     = false);

    //! Add an integrand which is just a single var. The a_var argument should
    //! correspond to the order in which the desired var was added to this
    //! object with add_var
    //! The last argument is whether to broadcast the result to all MPI ranks
    //! or just keep on rank 0. Most use cases won't need this set to true.
    void add_var_integrand(
        int a_var, std::vector<amrex::ParticleReal> &out_integrals,
        const IntegrationMethod &a_method_u = IntegrationMethod::trapezium,
        const IntegrationMethod &a_method_v = IntegrationMethod::trapezium,
        const bool a_broadcast_integral     = false);

    //! Integrate the integrands added using add_integrand
    void integrate();

    //! This integrate function can be used if you only want to integrate one
    //! integrand. It calls add_integrand() and integrate()
    //! The last argument is whether to broadcast the result to all MPI ranks
    //! or just keep on rank 0. Most use cases won't need this set to true.
    std::vector<amrex::ParticleReal> integrate(
        integrand_t a_integrand,
        const IntegrationMethod &a_method_u = IntegrationMethod::trapezium,
        const IntegrationMethod &a_method_v = IntegrationMethod::trapezium,
        const bool a_broadcast_integral     = false);

    //! Write the interpolated data to a file with a block for each surface
    void write_extraction(std::string a_file_prefix) const;

    //! write some integrals to a file at this timestep
    void write_integrals(
        const std::string &a_filename,
        const std::vector<std::vector<amrex::ParticleReal>> &a_integrals,
        const std::vector<std::string> &a_labels = {}) const;

    //! convenience caller for write_integrals in the case of just integral per
    //! surface
    void write_integral(const std::string &a_filename,
                        const std::vector<amrex::ParticleReal> &a_integrals,
                        const std::string &a_label = "") const;
};

#include "SurfaceExtraction.impl.hpp"

#endif /* SURFACEEXTRACTION_HPP_ */
