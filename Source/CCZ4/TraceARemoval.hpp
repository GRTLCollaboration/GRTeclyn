/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// This class enforces A to be trace-free
#ifndef TRACEAREMOVAL_HPP_
#define TRACEAREMOVAL_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"

class TraceARemoval
{
  public:
    struct Vars
    {
        Tensor<2, amrex::Real> h;
        Tensor<2, amrex::Real> A;

        template <typename mapping_function_t>
        AMREX_GPU_HOST_DEVICE void
        enum_mapping(mapping_function_t mapping_function);
    };

    AMREX_GPU_HOST_DEVICE void
    operator()(const amrex::CellData<amrex::Real> &cell) const
    {
        // NOLINTNEXTLINE [cppcoreguidelines-pro-type-member-init]
        Vars vars;
        load_vars(cell, vars);

        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
        TensorAlgebra::make_trace_free(vars.A, vars.h, h_UU);

        store_vars(cell, vars);
    }
};

template <typename mapping_function_t>
AMREX_GPU_HOST_DEVICE void
TraceARemoval::Vars::enum_mapping(mapping_function_t mapping_function)
{
    VarsTools::define_symmetric_enum_mapping(mapping_function,
                                             GRInterval<c_h11, c_h33>(), h);
    VarsTools::define_symmetric_enum_mapping(mapping_function,
                                             GRInterval<c_A11, c_A33>(), A);
}

#endif /* TRACEAREMOVAL_HPP_ */
