/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TRACEAREMOVAL_HPP_
#define TRACEAREMOVAL_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

// This class enforces det(h) = 1 and A to be trace-free
class TraceARemoval
{
  public:

    // Constructor
    TraceARemoval() = default;

    // Compute function
    AMREX_GPU_DEVICE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const
    {
        const amrex::CellData<amrex::Real> &state_cell_data =
            state.cellData(ix, iy, iz);
        // Copy here, to ensure reference in constructor stays in scope
        const amrex::CellData<const amrex::Real> &const_state_cell_data =
            state_cell_data;
        const CCZ4Vars vars(const_state_cell_data);

        using namespace CCZ4Geometry;

        // Enforce the unit determinant constraint on the conformal metric.
        const amrex::Real det_h = compute_metric_determinant(vars);
        AMREX_ASSERT(det_h > 0.0);
        const amrex::Real metric_factor =
            std::pow(det_h, -1.0 / static_cast<double>(GR_SPACEDIM));
        FOR2_SYM(i, j)
        {
            state_cell_data[VAR_IDX(c_h11, i, j)] *= metric_factor;
        }

        // vars references state_cell_data, so compute the trace using the
        // normalized conformal metric.
        const auto trace_A = compute_trace_A(vars);
        const amrex::Real one_over_gr_spacedim =
            1.0 / static_cast<double>(GR_SPACEDIM);
        FOR2_SYM(i, j)
        {
            state_cell_data[VAR_IDX(c_A11, i, j)] -=
                one_over_gr_spacedim * vars.h(i, j) * trace_A;
        }
    }
};

#endif /* TRACEAREMOVAL_HPP_ */
