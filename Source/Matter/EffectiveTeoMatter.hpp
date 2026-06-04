/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef EFFECTIVETEOMATTER_HPP_
#define EFFECTIVETEOMATTER_HPP_

#include "CCZ4AdvecVars.hpp"
#include "CCZ4D1Vars.hpp"
#include "CCZ4D2Vars.hpp"
#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "GRParmParse.hpp"
#include "TensorAlgebra.hpp"

//! State accessors for prescribed Teo effective stress-energy fields.
class EffectiveTeoMatterVars : public CCZ4Vars
{
  public:
    AMREX_GPU_DEVICE
    EffectiveTeoMatterVars(
        const amrex::CellData<const amrex::Real> &input_cell_data)
        : CCZ4Vars(input_cell_data)
    {
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    teo_rho() const
    {
        return cell_data[c_teo_rho];
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    teo_j(int i) const
    {
        return cell_data[c_teo_j1 + i];
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    teo_S(int i, int j) const
    {
        return cell_data[VAR_IDX(c_teo_S11, i, j)];
    }
};

//! Prescribed effective matter source implied by the Teo geometry.
class EffectiveTeoMatter
{
  protected:
    double m_source_strength;

  public:
    EffectiveTeoMatter() : m_source_strength(1.0)
    {
        GRParmParse pp;
        pp.load("teo_source_strength", m_source_strength, 1.0);
    }

    using Vars      = EffectiveTeoMatterVars;
    using D1Vars    = CCZ4D1Vars;
    using D2Vars    = CCZ4D2Vars;
    using AdvecVars = CCZ4AdvecVars;

    [[nodiscard]] AMREX_GPU_DEVICE emtensor_t
    compute_emtensor(const Vars &vars, const D1Vars &,
                     const Tensor<2, amrex::Real> &h_UU,
                     const Tensor<3, amrex::Real> &) const
    {
        emtensor_t out;
        out.rho = m_source_strength * vars.teo_rho();
        FOR (i)
        {
            out.j[i] = m_source_strength * vars.teo_j(i);
            FOR (j) { out.S[i][j] = m_source_strength * vars.teo_S(i, j); }
        }
        out.trS = vars.chi() * TensorAlgebra::compute_trace(out.S, h_UU);
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
        // The first Teo implementation prescribes an effective source; it does
        // not evolve a fundamental matter model.
    }
};

#endif /* EFFECTIVETEOMATTER_HPP_ */
