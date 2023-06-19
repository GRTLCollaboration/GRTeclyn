/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DERIVATIVETESTSCOMPUTE_HPP_
#define DERIVATIVETESTSCOMPUTE_HPP_

#include "Cell.hpp"
#include "VarsTools.hpp"

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

    template <class data_t> struct Vars
    {
        data_t d1;
        data_t d2;
        data_t d2_mixed;
        data_t diss;
        data_t advec_up;
        data_t advec_down;

        template <typename mapping_function_t>
        AMREX_GPU_DEVICE void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_d1, d1);
            define_enum_mapping(mapping_function, c_d2, d2);
            define_enum_mapping(mapping_function, c_d2_mixed, d2_mixed);
            define_enum_mapping(mapping_function, c_diss, diss);
            define_enum_mapping(mapping_function, c_advec_up, advec_up);
            define_enum_mapping(mapping_function, c_advec_down, advec_down);
        }
    };

    template <class data_t>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    operator()(int i, int j, int k, const amrex::Array4<data_t> &out,
               const amrex::Array4<data_t const> &in) const
    {
        const auto out_d1 = m_deriv.template diff1<Vars>(i, j, k, in);
        const auto out_d2 = m_deriv.template diff2<Vars>(i, j, k, in);

        Vars<data_t> out_diss;
        VarsTools::assign(out_diss, 0.);
        m_deriv.add_dissipation(i, j, k, out_diss, in, 1.0);

        Tensor<1, data_t> shift_down = {-2., 0., -3.};
        const auto out_advec_down =
            m_deriv.template advection<Vars>(i, j, k, in, shift_down);

        Tensor<1, data_t> shift_up = {2., 0., 3.};
        const auto out_advec_up =
            m_deriv.template advection<Vars>(i, j, k, in, shift_up);

        const auto out_cell_data = out.cellData(i, j, k);

        out_cell_data[c_d1]         = out_d1.d1[2];
        out_cell_data[c_d2]         = out_d2.d2[2][2];
        out_cell_data[c_d2_mixed]   = out_d2.d2[0][2];
        out_cell_data[c_diss]       = out_diss.diss;
        out_cell_data[c_advec_down] = out_advec_down.advec_down;
        out_cell_data[c_advec_up]   = out_advec_up.advec_up;
    }
};

#endif /* DERIVATIVETESTSCOMPUTE_HPP_ */
