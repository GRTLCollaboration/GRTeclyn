/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef NOMATTER_HPP_
#define NOMATTER_HPP_

#include "CCZ4AdvecVars.hpp"
#include "CCZ4D1Vars.hpp"
#include "CCZ4D2Vars.hpp"
#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"

//! Matter model for geometry-only controls.
class NoMatter
{
  public:
    using Vars      = CCZ4Vars;
    using D1Vars    = CCZ4D1Vars;
    using D2Vars    = CCZ4D2Vars;
    using AdvecVars = CCZ4AdvecVars;

    [[nodiscard]] AMREX_GPU_DEVICE emtensor_t
    compute_emtensor(const Vars &, const D1Vars &,
                     const Tensor<2, amrex::Real> &,
                     const Tensor<3, amrex::Real> &) const
    {
        emtensor_t out;
        out.rho  = 0.0;
        out.trS  = 0.0;
        FOR (i)
        {
            out.j[i] = 0.0;
            FOR (j) { out.S[i][j] = 0.0; }
        }
        return out;
    }

    [[nodiscard]] AMREX_GPU_DEVICE emtensor_t
    compute_emtensor(const Vars &vars, const D1Vars &d1,
                     const Tensor<2, amrex::Real> &h_UU,
                     const Tensor<3, amrex::Real> &chris_ULL,
                     const Coordinates &, amrex::Real) const
    {
        return compute_emtensor(vars, d1, h_UU, chris_ULL);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_matter_rhs(const amrex::CellData<amrex::Real> &, const Vars &,
                   const D1Vars &, const D2Vars &, const AdvecVars &) const
    {
    }
};

#endif /* NOMATTER_HPP_ */
