/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef DUSTMATTER_HPP_
#define DUSTMATTER_HPP_

#include "CCZ4AdvecVars.hpp"
#include "CCZ4D1Vars.hpp"
#include "CCZ4D2Vars.hpp"
#include "CCZ4Geometry.hpp"
#include "DustMatterVars.hpp"
#include "TensorAlgebra.hpp"

/// Pressureless dust with rest-mass density and 3-velocity (simplified advection).
class DustMatter
{
  public:
    using Vars = DustMatterVars;
    using D1Vars = CCZ4D1Vars;
    using D2Vars = CCZ4D2Vars;
    using AdvecVars = CCZ4AdvecVars;

    [[nodiscard]] AMREX_GPU_DEVICE emtensor_t
    compute_emtensor(const Vars &vars, const D1Vars &,
                     const Tensor<2, amrex::Real> &h_UU,
                     const Tensor<3, amrex::Real> &) const
    {
        emtensor_t out;
        const amrex::Real rho = amrex::max(vars.dust_rho(), amrex::Real(0.0));
        const amrex::Real alpha = vars.lapse();
        const amrex::Real inv_alpha = 1.0 / alpha;
        Tensor<1, amrex::Real> u_low;
        FOR (i)
        {
            u_low[i] = alpha * vars.dust_v(i) - vars.shift(i);
        }
        const amrex::Real u0 = -alpha;
        FOR (i, j)
        {
            out.S[i][j] = rho * u_low[i] * u_low[j] / (alpha * alpha);
        }
        out.rho = rho;
        out.trS = vars.chi() * TensorAlgebra::compute_trace(out.S, h_UU);
        FOR (i)
        {
            out.j[i] = rho * u0 * (vars.dust_v(i) + vars.shift(i) * inv_alpha);
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
    add_matter_rhs(const amrex::CellData<amrex::Real> &rhs, const Vars &,
                   const D1Vars &, const D2Vars &, const AdvecVars &) const
    {
        // Prescribed/co-moving dust for transport studies; full Valencia
        // conservative hydro evolution is future work.
        amrex::ignore_unused(rhs);
    }
};

#endif /* DUSTMATTER_HPP_ */
