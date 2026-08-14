/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef KLEINGORDONRHS_HPP_
#define KLEINGORDONRHS_HPP_

// GRTeclyn includes
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRParmParse.hpp"
#include "TensorAlgebra.hpp"

// Problem specific includes
#include "StateVariables.hpp"

template <class model_t, class deriv_t = FourthOrderDerivatives>
class KleinGordonRHS
{
  public:

    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    KleinGordonRHS(amrex::Real a_dx, model_t a_model)
        : m_deriv(a_dx), m_model(a_model)
    {
        GRParmParse pp;
        pp.get("grteclyn.sigma", m_sigma);
    };

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &rhs_state,
               const amrex::Array4<amrex::Real const> &state) const;

  private:
    amrex::Real m_sigma{};
    deriv_t m_deriv;
    model_t m_model;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void rhs_equation(
        const amrex::CellData<amrex::Real> &rhs_cell_data,
        const amrex::CellData<amrex::Real const> &state_cell_data,
        const amrex::Array1D<amrex::Real, 0, AMREX_SPACEDIM> &d2_phi) const;
};

#include "KleinGordonRHS.impl.hpp"

#endif
