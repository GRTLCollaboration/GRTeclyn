/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TRACEAREMOVAL_HPP_
#define TRACEAREMOVAL_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

// This class enforces A to be trace-free
class TraceARemoval
{
  public:

    // Constructor
    TraceARemoval() = default;

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
        const auto trace_A = compute_trace_A(vars);
        const amrex::Real one_over_gr_spacedim =
            1. / ((amrex::Real)GR_SPACEDIM);
        FOR2_SYM(i, j)
        {
            state_cell_data[sym_var_idx(c_A11, i, j)] -=
                one_over_gr_spacedim * vars.h(i, j) * trace_A;
        }
    }
};

#endif /* TRACEAREMOVAL_HPP_ */
