/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TRACEAREMOVAL_HPP_
#define TRACEAREMOVAL_HPP_

#include "CCZ4Vars2.hpp"
#include "Cell.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra2.hpp"

// This class enforces A to be trace-free
class TraceARemoval
{
  public:

    // Constructor
    TraceARemoval() {}

    // Compute function
    AMREX_GPU_HOST_DEVICE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const
    {
        const amrex::CellData<amrex::Real> &state_cell_data =
            state.cellData(ix, iy, iz);
        const CCZ4Vars2 vars(state_cell_data);

        using namespace TensorAlgebra2;
        const auto h_UU                   = compute_inverse_metric(vars);
        const auto trace_A                = compute_trace_A(vars, h_UU);
        const double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
        FOR (i, j)
        {
            state_cell_data[var_idx(c_A11, i, j)] -=
                one_over_gr_spacedim * vars.h(i, j) * trace_A;
        }
    }
};

#endif /* TRACEAREMOVAL_HPP_ */
