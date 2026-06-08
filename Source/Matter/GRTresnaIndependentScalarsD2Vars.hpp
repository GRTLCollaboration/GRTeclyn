#ifndef GRTRESNA_INDEPENDENT_SCALARS_D2VARS_HPP_
#define GRTRESNA_INDEPENDENT_SCALARS_D2VARS_HPP_

#include "CCZ4D2Vars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRTresnaScalarLayout.hpp"
#include "Tensor.hpp"

class GRTresnaIndependentScalarsD2Vars : public CCZ4D2Vars
{
  public:
    AMREX_GPU_DEVICE
    GRTresnaIndependentScalarsD2Vars(
        int ix, int iy, int iz, const amrex::Array4<const amrex::Real> &state,
        const FourthOrderDerivatives &a_deriv)
        : CCZ4D2Vars(ix, iy, iz, state, a_deriv)
    {
        for (int k = 0; k < GRTRESNA_MAX_INDEPENDENT_SCALARS; ++k)
        {
            m_scalar_d2_phi[k] =
                a_deriv.diff2_scalar(ix, iy, iz, state, c_phi_lump_index(k));
        }
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const Tensor<2, amrex::Real> &
    phi(int k) const
    {
        return m_scalar_d2_phi[k];
    }

  private:
    amrex::GpuArray<Tensor<2, amrex::Real>, GRTRESNA_MAX_INDEPENDENT_SCALARS>
        m_scalar_d2_phi{};
};

#endif /* GRTRESNA_INDEPENDENT_SCALARS_D2VARS_HPP_ */
