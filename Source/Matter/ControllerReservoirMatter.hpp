#ifndef CONTROLLER_RESERVOIR_MATTER_HPP_
#define CONTROLLER_RESERVOIR_MATTER_HPP_

#include "CCZ4Geometry.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "RLMatterPumpParams.hpp"
#include "RLPumpForce.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

//! Decorator: wraps an InnerMatter and adds a controller energy-momentum
//! reservoir (ρ_c, j_{c i}) that absorbs the PD-pump 4-force.
//!
//! Modes (controller_reservoir_mode):
//!   0 — inactive (bit-identical to InnerMatter alone)
//!   1 — ledger: evolve reservoir; include in compute_emtensor when
//!       include_emtensor=true (constraints path only)
//!   2 — backreaction: evolve reservoir; include in compute_emtensor whenever
//!       include_emtensor=true (constraints + RHS paths)
//!
//! IsBicomplex selects the field layout / force-density accessors at
//! compile time so BiComplex and single-complex D1Vars stay incompatible.
template <class InnerMatter, bool IsBicomplex>
class ControllerReservoirMatter
{
  public:
    using Vars      = typename InnerMatter::Vars;
    using D1Vars    = typename InnerMatter::D1Vars;
    using D2Vars    = typename InnerMatter::D2Vars;
    using AdvecVars = typename InnerMatter::AdvecVars;

    ControllerReservoirMatter(InnerMatter a_inner, RLMatterPumpParams a_pump,
                              int a_mode, bool a_include_emtensor)
        : m_inner(a_inner), m_pump(a_pump), m_mode(a_mode),
          m_include_emtensor(a_include_emtensor)
    {
    }

    [[nodiscard]] AMREX_GPU_DEVICE emtensor_t compute_emtensor(
        const Vars &vars, const D1Vars &d1, const Tensor<2, amrex::Real> &h_UU,
        const Tensor<3, amrex::Real> &chris_ULL) const
    {
        auto em = m_inner.compute_emtensor(vars, d1, h_UU, chris_ULL);
        add_reservoir_to_emtensor(em, vars);
        return em;
    }

    [[nodiscard]] AMREX_GPU_DEVICE emtensor_t
    compute_emtensor(const Vars &vars, const D1Vars &d1,
                     const Tensor<2, amrex::Real> &h_UU,
                     const Tensor<3, amrex::Real> &chris_ULL,
                     const Coordinates &coords, amrex::Real time) const
    {
        auto em =
            m_inner.compute_emtensor(vars, d1, h_UU, chris_ULL, coords, time);
        add_reservoir_to_emtensor(em, vars);
        return em;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_matter_rhs(const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
                   const D1Vars &d1, const D2Vars &d2,
                   const AdvecVars &advec) const
    {
        m_inner.add_matter_rhs(rhs, vars, d1, d2, advec);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_matter_rhs(const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
                   const D1Vars &d1, const D2Vars &d2, const AdvecVars &advec,
                   const Coordinates &coords, amrex::Real time) const
    {
        m_inner.add_matter_rhs(rhs, vars, d1, d2, advec, coords, time);
        if (m_mode < 1)
        {
            return;
        }
        add_reservoir_rhs(rhs, vars, d1, coords, time);
    }

  private:
    InnerMatter m_inner;
    RLMatterPumpParams m_pump{};
    int m_mode{0};
    bool m_include_emtensor{false};

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_reservoir_to_emtensor(emtensor_t &em, const Vars &vars) const
    {
        if (!m_include_emtensor || m_mode < 1)
        {
            return;
        }
        em.rho += vars.cell_data[c_rho_ctrl];
        em.j[0] += vars.cell_data[c_jctrl1];
        em.j[1] += vars.cell_data[c_jctrl2];
        em.j[2] += vars.cell_data[c_jctrl3];
    }

    //! Local reservoir evolution absorbing the pump 4-force:
    //!   ∂_t ρ_c = α K ρ_c − α f_⊥
    //!   ∂_t j_i = −ρ_c ∂_i α − α f_i
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_reservoir_rhs(const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
                      const D1Vars &d1, const Coordinates &coords,
                      amrex::Real time) const
    {
        const amrex::Real lapse = vars.lapse();
        const amrex::Real K     = vars.K();
        const amrex::Real rho_c = vars.cell_data[c_rho_ctrl];

        const RLPumpForceDensity force = compute_force(vars, d1, coords, time);

        rhs[c_rho_ctrl] = lapse * K * rho_c - lapse * force.f_perp;
        rhs[c_jctrl1]   = -rho_c * d1.lapse(0) - lapse * force.f_x;
        rhs[c_jctrl2]   = -rho_c * d1.lapse(1) - lapse * force.f_y;
        rhs[c_jctrl3]   = -rho_c * d1.lapse(2) - lapse * force.f_z;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE RLPumpForceDensity
    compute_force(const Vars &vars, const D1Vars &d1, const Coordinates &coords,
                  amrex::Real time) const
    {
        const amrex::Real lapse = vars.lapse();
        if constexpr (IsBicomplex)
        {
            const amrex::Real phi1p = vars.cell_data[c_phi];
            const amrex::Real Pi1p  = vars.cell_data[c_Pi];
            const amrex::Real phi2p = vars.cell_data[c_phi2];
            const amrex::Real Pi2p  = vars.cell_data[c_Pi2];
            const amrex::Real phi1m = vars.cell_data[c_phi_m];
            const amrex::Real Pi1m  = vars.cell_data[c_Pi_m];
            const amrex::Real phi2m = vars.cell_data[c_phi2_m];
            const amrex::Real Pi2m  = vars.cell_data[c_Pi2_m];
            const RLPumpSources src = RLPumpForce::compute_bicomplex_sources(
                m_pump, coords.x, coords.y, coords.z, time, lapse, phi1p, phi2p,
                Pi1p, Pi2p, phi1m, phi2m, Pi1m, Pi2m);
            return RLPumpForce::force_density_from_sources(
                src, 1.0, Pi1p, Pi2p, d1.phi1p(0), d1.phi1p(1), d1.phi1p(2),
                d1.phi2p(0), d1.phi2p(1), d1.phi2p(2), -1.0, Pi1m, Pi2m,
                d1.phi1m(0), d1.phi1m(1), d1.phi1m(2), d1.phi2m(0), d1.phi2m(1),
                d1.phi2m(2));
        }
        else
        {
            const amrex::Real phi1 = vars.cell_data[c_phi];
            const amrex::Real Pi1  = vars.cell_data[c_Pi];
            const amrex::Real phi2 = vars.cell_data[c_phi2];
            const amrex::Real Pi2  = vars.cell_data[c_Pi2];
            const RLPumpSources src = RLPumpForce::compute_single_field_sources(
                m_pump, coords.x, coords.y, coords.z, time, lapse, phi1, phi2,
                Pi1, Pi2);
            return RLPumpForce::force_density_single(
                src, 1.0, Pi1, Pi2, d1.phi1(0), d1.phi1(1), d1.phi1(2),
                d1.phi2(0), d1.phi2(1), d1.phi2(2));
        }
    }
};

#endif /* CONTROLLER_RESERVOIR_MATTER_HPP_ */
