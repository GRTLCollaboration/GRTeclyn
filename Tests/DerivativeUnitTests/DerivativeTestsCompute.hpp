/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef DERIVATIVETESTSCOMPUTE_HPP_
#define DERIVATIVETESTSCOMPUTE_HPP_

#include "Cell.hpp"
#include "Tensor.hpp"

enum
{
    c_d1,
    c_d2,
    c_d2_mixed,
    c_diss,
    c_advec_up,
    c_advec_down,
    NUM_DERIVATIVES_VARS
};

template <class deriv_t> class DerivativeTestsCompute
{
  private:
    deriv_t m_deriv;

  public:

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE DerivativeTestsCompute(double dx)
        : m_deriv(dx)
    {
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz, const amrex::Array4<amrex::Real> &out,
               const amrex::Array4<amrex::Real const> &in) const
    {
        TensorArray::Rank1 out_d1 =
            m_deriv.diff1_array_scalar(ix, iy, iz, in, c_d1);
        TensorArray::Rank2 out_d2 =
            m_deriv.diff2_array_scalar(ix, iy, iz, in, c_d2);

        TensorArray::Rank1 shift_up   = {2., 0., 3.};
        TensorArray::Rank1 shift_down = {-2., 0., -3.};

        const amrex::Real out_advec_down =
            m_deriv.advection(ix, iy, iz, in, shift_down, c_advec_down);

        const amrex::Real out_advec_up =
            m_deriv.advection(ix, iy, iz, in, shift_up, c_advec_up);

        double sigma = 1.0;
        amrex::Real out_diss =
            m_deriv.calculate_dissipation(ix, iy, iz, in, sigma, c_diss);

        const auto out_cell_data = out.cellData(ix, iy, iz);

        out_cell_data[c_d1]         = out_d1(2);
        out_cell_data[c_d2]         = out_d2(2, 2);
        out_cell_data[c_d2_mixed]   = out_d2(0, 2);
        out_cell_data[c_diss]       = out_diss;
        out_cell_data[c_advec_down] = out_advec_down;
        out_cell_data[c_advec_up]   = out_advec_up;
    }
};

#endif /* DERIVATIVETESTSCOMPUTE_HPP_ */
