/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef ALGEBRAICCONSTRAINTSENFORCER_HPP_
#define ALGEBRAICCONSTRAINTSENFORCER_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

using namespace amrex::literals;

// This class enforces det(h)=1 and A to be trace-free
class AlgebraicConstraintsEnforcer
{
  public:

    // Constructor
    AlgebraicConstraintsEnforcer() = default;

    // Compute function

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
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
        AMREX_ASSERT(det_h > 0.0_rt);
        const amrex::Real metric_factor =
            std::pow(det_h, -1.0_rt / static_cast<amrex::Real>(GR_SPACEDIM));
        FOR2_SYM(i, j)
        {
            state_cell_data[sym_var_idx(c_h11, i, j)] *= metric_factor;
        }
        // Enforce A to be trace-free
        const auto trace_A = compute_trace_A(vars);
        const amrex::Real one_over_gr_spacedim =
            1._rt / ((amrex::Real)GR_SPACEDIM);
        FOR2_SYM(i, j)
        {
            state_cell_data[sym_var_idx(c_A11, i, j)] -=
                one_over_gr_spacedim * vars.h(i, j) * trace_A;
        }
    }
};

#endif /* ALGEBRAICCONSTRAINTSENFORCER_HPP_ */
