/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4RHS_HPP_
#define CCZ4RHS_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MovingPunctureGauge.hpp"
#include "TensorAlgebra.hpp"

#include "StateVariables.hpp" //This files needs NUM_VARS - total number of components

#include <array>

/// Base parameter struct for CCZ4
/** This struct collects the gauge independent CCZ4 parameters i.e. the damping
 * ones
 */
struct CCZ4_base_params_t
{
    double kappa1;    //!< Damping parameter kappa1 as in arXiv:1106.2254
    double kappa2;    //!< Damping parameter kappa2 as in arXiv:1106.2254
    double kappa3;    //!< Damping parameter kappa3 as in arXiv:1106.2254
    bool covariantZ4; //!< if true, replace kappa1->kappa1/lapse as in
                      //!<  arXiv:1307.7391 eq. 27
};

/// Parameter struct for CCZ4
/** This struct collects all parameters that are necessary for CCZ4 such as
 * gauge and damping parameters. It inherits from CCZ4_base_params_t and
 * gauge_t::params_t
 */
template <class gauge_params_t = MovingPunctureGauge::params_t>
struct CCZ4_params_t : public CCZ4_base_params_t, public gauge_params_t
{
};

/// Compute class to calculate the CCZ4 right hand side
/**
 * This compute class implements the CCZ4 right hand side equations.
 **/
template <class gauge_t = MovingPunctureGauge,
          class deriv_t = FourthOrderDerivatives>
class CCZ4RHS
{
  public:
    enum formulations : int
    {
        USE_CCZ4 = 0,
        USE_BSSN = 1
    };

    enum covariantZ4 : int
    {
        YES,
        NO
    };

    using params_t = CCZ4_params_t<typename gauge_t::params_t>;

  protected:
    params_t m_params; //!< CCZ4 parameters
    gauge_t m_gauge;   //!< Class to compute gauge in rhs_equation
    double m_sigma;    //!< Coefficient for Kreiss-Oliger dissipation
    int m_formulation;
    double m_cosmological_constant;
    deriv_t m_deriv;

  public:
    /// Constructor
    CCZ4RHS(
        params_t a_params,            //!< The CCZ4 parameters
        double a_dx,                  //!< The grid spacing
        double a_sigma,               //!< Kreiss-Oliger dissipation coefficient
        int a_formulation = USE_CCZ4, //!< Switches between CCZ4, BSSN,...
        double a_cosmological_constant = 0 //!< Value of the cosmological const.
    );

    /// Compute function
    /** This function orchestrates the calculation of the rhs for one specific
     * grid cell.
     */
    /// Calculates the rhs for chi and h_ij
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_chi_and_h_ij(int ix, int iy, int iz,
                         const amrex::Array4<amrex::Real> &rhs,
                         const amrex::Array4<const amrex::Real> &state) const;

    // Calculates rhs for A_ij and Theta and Gamma
    template <int formulation, int use_covariant_Z4>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_A_ij_and_Theta_and_Gamma(
        int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs,
        const amrex::Array4<const amrex::Real> &state) const;

    // Apply gauage (no derivatives needed here!)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    calculate_gauge_rhs(int ix, int iy, int iz,
                        const amrex::Array4<amrex::Real> &rhs,
                        const amrex::Array4<const amrex::Real> &state) const;

    // Apply dissipation (split for matter classes)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    apply_dissipation(int ix, int iy, int iz,
                      const amrex::Array4<amrex::Real> &rhs,
                      const amrex::Array4<const amrex::Real> &state) const;
};

#include "CCZ4RHS.impl.hpp"

#endif /* CCZ4RHS_HPP_ */
