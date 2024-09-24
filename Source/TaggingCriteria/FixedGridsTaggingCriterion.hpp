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
    double m_dx;
    double m_L;
    int m_level;
    std::array<double, AMREX_SPACEDIM> m_center;

  public:
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    FixedGridsTaggingCriterion(
        const double dx, const int a_level, const double a_L,
        const std::array<double, AMREX_SPACEDIM> a_center)
        : m_dx(dx), m_L(a_L), m_level(a_level), m_center(a_center){};
    // NOLINTEND(bugprone-easily-swappable-parameters)
    template <class data_t>
    AMREX_GPU_DEVICE data_t
    compute(int i, int j, int k,
            const amrex::Array4<const data_t> &current_cell) const
    {
        data_t criterion = 0.0;
        // make sure the inner part is regridded around the horizon
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4
        double ratio = pow(2.0, -(m_level + 2.0));

        amrex::IntVect cell(AMREX_D_DECL(i, j, k));

        const Coordinates<data_t> coords(cell, m_dx, m_center);
        const data_t max_abs_xy =
            simd_max(std::abs(coords.x), std::abs(coords.y));
        const data_t max_abs_xyz = simd_max(max_abs_xy, std::abs(coords.z));
        auto regrid              = simd_compare_lt(max_abs_xyz, m_L * ratio);
        criterion                = simd_conditional(regrid, 100.0, criterion);

        return criterion;
    }
};

#endif /* FIXEDGRIDSTAGGINGCRITERION_HPP_ */
