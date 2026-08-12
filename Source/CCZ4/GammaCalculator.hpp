/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef GAMMACALCULATOR_HPP_
#define GAMMACALCULATOR_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"

class GammaCalculator
{
  public:
    AMREX_GPU_HOST_DEVICE
        AMREX_FORCE_INLINE explicit GammaCalculator(double a_dx)
        : m_deriv(a_dx)
    {
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const
    {
        const amrex::CellData<amrex::Real> &state_cell_data =
            state.cellData(ix, iy, iz);
        const amrex::CellData<const amrex::Real> &const_state_cell_data =
            state_cell_data;
        const CCZ4Vars vars(const_state_cell_data);

        const auto h_UU = CCZ4Geometry::compute_inverse_metric(vars);
        const auto d1_h = m_deriv.d1_sym_tensor(ix, iy, iz, state, c_h11);
        const auto christoffel = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

        FOR (i)
        {
            state_cell_data[c_Gamma1 + i] = christoffel.contracted(i);
        }
    }

  private:
    FourthOrderDerivatives m_deriv;
};

#endif /* GAMMACALCULATOR_HPP_ */
