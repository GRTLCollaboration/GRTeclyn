/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef KLEINGORDONRHS_HPP_
#define KLEINGORDONRHS_HPP_

// GRTeclyn includes
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"

// Problem specific includes
#include "StateVariables.hpp"

template <class model_t, class deriv_t = FourthOrderDerivatives>
class KleinGordonRHS
{
  public:

    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    KleinGordonRHS(amrex::Real a_sigma, amrex::Real a_dx, model_t a_model)
        : m_sigma(a_sigma), m_deriv(a_dx), m_model(a_model) {};

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute(int i, int j, int k, const amrex::Array4<data_t const> &input,
            const amrex::Array4<data_t> &output) const;

  private:
    amrex::Real m_sigma;
    deriv_t m_deriv;
    model_t m_model;

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    rhs_equation(const amrex::CellData<data_t const> &input_cell_data,
                 const amrex::CellData<data_t> &output_cell_data,
                 const amrex::Array1D<data_t, 0, AMREX_SPACEDIM> &d2phi) const;
};

#include "KleinGordonRHS.impl.hpp"

#endif
