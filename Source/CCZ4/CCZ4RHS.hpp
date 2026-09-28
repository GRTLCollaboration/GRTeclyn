/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4RHS_HPP_
#define CCZ4RHS_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRParmParse.hpp"
#include "MovingPunctureGauge.hpp"
#include "TensorAlgebra.hpp"

#include "StateVariables.hpp" //This files needs NUM_VARS - total number of components

#include <array>

using namespace amrex::literals;

/// Base parameter struct for CCZ4
/** This struct collects the gauge-independent CCZ4 parameters. */
struct CCZ4_params_t
{
    amrex::Real bssn_coeff; //!< 0 for CCZ4 and 1 for BSSN
    amrex::Real kappa1;     //!< Damping parameter kappa1 as in arXiv:1106.2254
    amrex::Real kappa2;     //!< Damping parameter kappa2 as in arXiv:1106.2254
    amrex::Real kappa3;     //!< Damping parameter kappa3 as in arXiv:1106.2254
    amrex::Real covariant_z4_coeff; //!< 1 to replace kappa1->kappa1/lapse as in
                                    //!< arXiv:1307.7391 eq. 27

    static void check_params();

    void fill_params()
    {
        GRParmParse ccz4_pp("ccz4");
        int formulation{};
        ccz4_pp.get("formulation", formulation);
        bssn_coeff = static_cast<amrex::Real>(formulation);
        ccz4_pp.get("kappa1", kappa1);
        ccz4_pp.get("kappa2", kappa2);
        ccz4_pp.get("kappa3", kappa3);
        bool covariant_z4{};
        ccz4_pp.get("covariantZ4", covariant_z4);
        covariant_z4_coeff = static_cast<amrex::Real>(covariant_z4);
    }
};

/// Compute class to calculate the CCZ4 right hand side
/**
 * This compute class implements the CCZ4 right hand side equations.
 **/
template <class deriv_t = FourthOrderDerivatives> class CCZ4RHS
{
  public:
    enum formulations : int
    {
        USE_CCZ4 = 0,
        USE_BSSN = 1
    };

    using params_t = CCZ4_params_t;

  protected:
    params_t m_params;   //!< CCZ4 parameters
    amrex::Real m_sigma; //!< Coefficient for Kreiss-Oliger dissipation
    amrex::Real m_cosmological_constant;
    deriv_t m_deriv;

  public:
    /// Constructor
    CCZ4RHS(amrex::Real a_dx, //!< The grid spacing
            amrex::Real a_cosmological_constant =
                0 //!< Value of the cosmological const.
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
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_A_ij_and_Theta_and_Gamma(
        int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs,
        const amrex::Array4<const amrex::Real> &state) const;

    // Apply dissipation (split for matter classes)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    apply_dissipation(int ix, int iy, int iz,
                      const amrex::Array4<amrex::Real> &rhs,
                      const amrex::Array4<const amrex::Real> &state) const;
};

inline void CCZ4_params_t::check_params()
{
    GRParmParse ccz4_pp("ccz4");

    int formulation{};
    formulation = CCZ4RHS<>::USE_CCZ4;
    ccz4_pp.queryAdd("formulation", formulation);
    if (formulation != CCZ4RHS<>::USE_CCZ4 &&
        formulation != CCZ4RHS<>::USE_BSSN)
    {
        ccz4_pp.error("formulation", "must be 0 (CCZ4) or 1 (BSSN)");
    }

    if (formulation == CCZ4RHS<>::USE_BSSN)
    {
        for (const char *kappa_name : {"kappa1", "kappa2", "kappa3"})
        {
            if (ccz4_pp.contains(kappa_name))
            {
                ccz4_pp.warning(
                    kappa_name,
                    "should not be provided with BSSN formulation, setting "
                    "it to zero");
            }
        }
        ccz4_pp.add("kappa1", 0.0);
        ccz4_pp.add("kappa2", 0.0);
        ccz4_pp.add("kappa3", 0.0);
    }
    else
    {
        amrex::Real kappa1 = 0.1_rt;
        ccz4_pp.queryAdd("kappa1", kappa1);
        if (kappa1 <= 0.0_rt)
        {
            ccz4_pp.warning("kappa1",
                            "should be greater than 0.0 to damp constraints "
                            "(see arXiv:1106.2254).");
        }

        amrex::Real kappa2 = 0.0_rt;
        ccz4_pp.queryAdd("kappa2", kappa2);
        if (kappa2 <= -1.0_rt)
        {
            ccz4_pp.warning("kappa2",
                            "should be greater than -1.0 to damp constraints "
                            "(see arXiv:1106.2254).");
        }

        amrex::Real kappa3 = 1.0_rt;
        ccz4_pp.queryAdd("kappa3", kappa3);
    }

    bool covariantZ4 = true;
    ccz4_pp.queryAdd("covariantZ4", covariantZ4);
}

#include "CCZ4RHS.impl.hpp"

#endif /* CCZ4RHS_HPP_ */
