/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef FIXEDGRIDSTAGGINGCRITERION_HPP_
#define FIXEDGRIDSTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

class FixedGridsTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const int m_level;
    const std::array<double, AMREX_SPACEDIM> m_center;

  public:
    FixedGridsTaggingCriterion(
        const amrex::Real dx, const int a_level, const amrex::Real a_L,
        const std::array<double, AMREX_SPACEDIM> a_center)
        : m_dx(dx), m_L(a_L), m_level(a_level), m_center(a_center){};

    template <class data_t>
    AMREX_GPU_DEVICE data_t compute(int i, int j, int k,
                                    const amrex::Array4<data_t> &state) const
    {
        data_t yes_regrid = 100.0;
        data_t no_regrid  = 0.0;
        // make sure the inner part is regridded around the horizon
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4
        double ratio = pow(2.0, -(m_level + 2.0));

        amrex::Real x = (i + 0.5) * m_dx - m_center[0];
        amrex::Real y = (j + 0.5) * m_dx - m_center[1];
        amrex::Real z = (k + 0.5) * m_dx - m_center[2];

        //	const Coordinates<data_t> coords(amrex::IntVect(AMREX_D_DECL(i,
        // j, k)), m_dx);
        const data_t max_abs_xy = simd_max(std::abs(x), std::abs(y));

        const data_t max_abs_xyz = simd_max(max_abs_xy, std::abs(z));
        auto regrid              = simd_compare_lt(max_abs_xyz, m_L * ratio);
        data_t criterion = simd_conditional(regrid, yes_regrid, no_regrid);

        // Write back into the flattened Chombo box
        //        store_vars(state.cellData(i,j,k), criterion);

        return criterion;
    }
};

#endif /* FIXEDGRIDSTAGGINGCRITERION_HPP_ */
