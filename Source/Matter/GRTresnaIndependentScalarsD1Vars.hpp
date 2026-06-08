#ifndef GRTRESNA_INDEPENDENT_SCALARS_D1VARS_HPP_
#define GRTRESNA_INDEPENDENT_SCALARS_D1VARS_HPP_

#include "CCZ4D1Vars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRTresnaIndependentScalarsVars.hpp"
#include "GRTresnaScalarLayout.hpp"
#include "Tensor.hpp"

class GRTresnaIndependentScalarsD1Vars : public CCZ4D1Vars
{
  public:
    AMREX_GPU_DEVICE
    GRTresnaIndependentScalarsD1Vars(
        int ix, int iy, int iz, const amrex::Array4<const amrex::Real> &state,
        const FourthOrderDerivatives &a_deriv)
        : CCZ4D1Vars(ix, iy, iz, state, a_deriv)
    {
        for (int k = 0; k < GRTRESNA_MAX_INDEPENDENT_SCALARS; ++k)
        {
            m_scalar_d1_phi[k] =
                a_deriv.diff1_scalar(ix, iy, iz, state, c_phi_lump_index(k));
            m_scalar_d1_Pi[k] =
                a_deriv.diff1_scalar(ix, iy, iz, state, c_Pi_lump_index(k));
        }
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &phi(int k,
                                                               int i) const
    {
        return m_scalar_d1_phi[k][i];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Pi(int k,
                                                                int i) const
    {
        return m_scalar_d1_Pi[k][i];
    }

  private:
    amrex::GpuArray<Tensor<1, amrex::Real>, GRTRESNA_MAX_INDEPENDENT_SCALARS>
        m_scalar_d1_phi{};
    amrex::GpuArray<Tensor<1, amrex::Real>, GRTRESNA_MAX_INDEPENDENT_SCALARS>
        m_scalar_d1_Pi{};
};

#endif /* GRTRESNA_INDEPENDENT_SCALARS_D1VARS_HPP_ */
