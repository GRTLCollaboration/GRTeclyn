#ifndef GRTRESNA_INDEPENDENT_SCALARS_ADVECVARS_HPP_
#define GRTRESNA_INDEPENDENT_SCALARS_ADVECVARS_HPP_

#include "CCZ4AdvecVars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRTresnaScalarLayout.hpp"
#include "Tensor.hpp"

class GRTresnaIndependentScalarsAdvecVars : public CCZ4AdvecVars
{
  public:
    AMREX_GPU_DEVICE
    GRTresnaIndependentScalarsAdvecVars(
        int ix, int iy, int iz, const amrex::Array4<const amrex::Real> &state,
        const FourthOrderDerivatives &a_deriv)
        : CCZ4AdvecVars(ix, iy, iz, state, a_deriv)
    {
        for (int k = 0; k < GRTRESNA_MAX_INDEPENDENT_SCALARS; ++k)
        {
            m_scalar_advec_phi[k] = a_deriv.advec_scalar(
                ix, iy, iz, state, m_shift_vector, c_phi_lump_index(k));
            m_scalar_advec_Pi[k] = a_deriv.advec_scalar(
                ix, iy, iz, state, m_shift_vector, c_Pi_lump_index(k));
        }
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &phi(int k) const
    {
        return m_scalar_advec_phi[k];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Pi(int k) const
    {
        return m_scalar_advec_Pi[k];
    }

  private:
    amrex::GpuArray<amrex::Real, GRTRESNA_MAX_INDEPENDENT_SCALARS>
        m_scalar_advec_phi{};
    amrex::GpuArray<amrex::Real, GRTRESNA_MAX_INDEPENDENT_SCALARS>
        m_scalar_advec_Pi{};
};

#endif /* GRTRESNA_INDEPENDENT_SCALARS_ADVECVARS_HPP_ */
