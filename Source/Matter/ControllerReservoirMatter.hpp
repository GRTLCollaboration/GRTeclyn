#ifndef CONTROLLER_RESERVOIR_MATTER_HPP_
#define CONTROLLER_RESERVOIR_MATTER_HPP_

#include "CCZ4Geometry.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "RLMatterPumpParams.hpp"
#include "RLPumpForce.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

//! Decorator: wraps an InnerMatter and adds a controller energy-momentum
//! reservoir (ρ_c, j_{c i}) that absorbs the PD-pump 4-force, so that
//! ∇_μ (T_matter + T_c)^{μν} = 0 up to truncation error.
//!
//! Modes (controller_reservoir_mode):
//!   0 — inactive (this decorator is not instantiated at all)
//!   1 — ledger: evolve reservoir; include in the constraint EMT only
//!   2 — backreaction: as 1, and also source the CCZ4 RHS from it
//! The mode-1 / mode-2 distinction is made by the CALLER, which passes
//! include_emtensor = (mode >= 1) on the constraints path and
//! include_emtensor = (mode >= 2) on the RHS path.
//!
//! Sign conventions (this codebase, verified against the scalar EMT):
//!   ρ      = ½ Σ_A s_A (Π_A² + |∇φ_A|²) + V     [BiComplexScalarField.impl]
//!   j_i    = Σ_A s_A (−Π_A ∂_i φ_A)
//! The pump adds S_A to ∂_t Π_A (a coordinate-time rate, no lapse), hence
//!   ∂_t ρ|_pump   = + Σ_A s_A Π_A S_A       ≡ +f_⊥
//!   ∂_t j_i|_pump = − Σ_A s_A S_A ∂_i φ_A   ≡ −f_i
//! Note the OPPOSITE relative sign: the reservoir must therefore absorb
//! −f_⊥ in the energy equation and +f_i in the momentum equation, and
//! neither source carries a lapse factor.
//!
//! Transport (S_c^{ij} = 0, i.e. a pressureless ledger):
//!   ∂_t ρ_c    = β^k∂_k ρ_c + αKρ_c − α D_i j_c^i − 2 j_c^i ∂_i α − f_⊥
//!   ∂_t j_{ci} = β^k∂_k j_{ci} + αK j_{ci} + j_{ck} ∂_i β^k − ρ_c ∂_i α + f_i
//! which is the standard 3+1 matter conservation system with the pump force
//! as source. Dropping the transport terms (as an earlier revision did) makes
//! ∇_μ T_c^{μν} = −f^ν false and the mode-2 Bianchi argument invalid.
//!
//! IsBicomplex selects the field layout / force-density accessors at
//! compile time so BiComplex and single-complex D1Vars stay incompatible.
template <class InnerMatter, bool IsBicomplex>
class ControllerReservoirMatter
{
  public:
    using Vars       = typename InnerMatter::Vars;
    using InnerD1    = typename InnerMatter::D1Vars;
    using InnerAdvec = typename InnerMatter::AdvecVars;
    using D2Vars     = typename InnerMatter::D2Vars;

    //! Inner D1Vars plus first derivatives of the four reservoir components.
    class D1Vars : public InnerD1
    {
      public:
        // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
        AMREX_GPU_DEVICE
        D1Vars(int ix, int iy, int iz,
               const amrex::Array4<const amrex::Real> &state,
               const FourthOrderDerivatives &a_deriv)
            : InnerD1(ix, iy, iz, state, a_deriv)
        {
            m_rho_ctrl_d1 =
                a_deriv.diff1_scalar(ix, iy, iz, state, c_rho_ctrl);
            const Tensor<1, amrex::Real> dj1 =
                a_deriv.diff1_scalar(ix, iy, iz, state, c_jctrl1);
            const Tensor<1, amrex::Real> dj2 =
                a_deriv.diff1_scalar(ix, iy, iz, state, c_jctrl2);
            const Tensor<1, amrex::Real> dj3 =
                a_deriv.diff1_scalar(ix, iy, iz, state, c_jctrl3);
            FOR (i)
            {
                m_jctrl_d1[0][i] = dj1[i];
                m_jctrl_d1[1][i] = dj2[i];
                m_jctrl_d1[2][i] = dj3[i];
            }
        }
        // NOLINTEND(cppcoreguidelines-pro-type-member-init)

        Tensor<1, amrex::Real> m_rho_ctrl_d1; //!< ∂_i ρ_c
        Tensor<2, amrex::Real> m_jctrl_d1;    //!< m_jctrl_d1[k][i] = ∂_i j_{ck}

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
        rho_ctrl(int i) const
        {
            return m_rho_ctrl_d1[i];
        }
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
        jctrl(int k, int i) const
        {
            return m_jctrl_d1[k][i];
        }
    };

    //! Inner AdvecVars plus upwinded β^k∂_k of the reservoir components.
    class AdvecVars : public InnerAdvec
    {
      public:
        // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
        AMREX_GPU_DEVICE
        AdvecVars(int ix, int iy, int iz,
                  const amrex::Array4<const amrex::Real> &state,
                  const FourthOrderDerivatives &a_deriv)
            : InnerAdvec(ix, iy, iz, state, a_deriv)
        {
            m_rho_ctrl_advec = a_deriv.advec_scalar(
                ix, iy, iz, state, this->m_shift_vector, c_rho_ctrl);
            m_jctrl_advec[0] = a_deriv.advec_scalar(
                ix, iy, iz, state, this->m_shift_vector, c_jctrl1);
            m_jctrl_advec[1] = a_deriv.advec_scalar(
                ix, iy, iz, state, this->m_shift_vector, c_jctrl2);
            m_jctrl_advec[2] = a_deriv.advec_scalar(
                ix, iy, iz, state, this->m_shift_vector, c_jctrl3);
        }
        // NOLINTEND(cppcoreguidelines-pro-type-member-init)

        amrex::Real m_rho_ctrl_advec;
        Tensor<1, amrex::Real> m_jctrl_advec;

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
        rho_ctrl() const
        {
            return m_rho_ctrl_advec;
        }
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
        jctrl(int k) const
        {
            return m_jctrl_advec[k];
        }
    };

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
        add_reservoir_rhs(rhs, vars, d1, advec, coords, time);
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
        // S_c^{ij} = 0, so only rho and j_i are contributed (trS untouched).
        em.rho += vars.cell_data[c_rho_ctrl];
        em.j[0] += vars.cell_data[c_jctrl1];
        em.j[1] += vars.cell_data[c_jctrl2];
        em.j[2] += vars.cell_data[c_jctrl3];
    }

    //! Reservoir evolution: standard 3+1 conservation with S_c^{ij} = 0 and
    //! the pump 4-force as source. See the class comment for the derivation
    //! and for why the momentum source is +f_i while the energy source is
    //! −f_⊥, and why neither carries a lapse.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_reservoir_rhs(const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
                      const D1Vars &d1, const AdvecVars &advec,
                      const Coordinates &coords, amrex::Real time) const
    {
        const amrex::Real lapse = vars.lapse();
        const amrex::Real K     = vars.K();
        const amrex::Real chi   = vars.chi();
        const amrex::Real rho_c = vars.cell_data[c_rho_ctrl];

        Tensor<1, amrex::Real> j_c;
        j_c[0] = vars.cell_data[c_jctrl1];
        j_c[1] = vars.cell_data[c_jctrl2];
        j_c[2] = vars.cell_data[c_jctrl3];

        const auto h_UU = CCZ4Geometry::compute_inverse_metric(vars);

        // Physical inverse spatial metric γ^{ij} = χ h^{ij}; raise the index.
        Tensor<1, amrex::Real> j_cU;
        FOR (i)
        {
            j_cU[i] = 0.0;
            FOR (j) { j_cU[i] += chi * h_UU[i][j] * j_c[j]; }
        }

        // D_i j_c^i = ∂_i j_c^i + j_c^i ∂_i ln√γ, with √γ = χ^{-3/2}
        // (det h = 1), so ∂_i ln√γ = −(3/2) ∂_i χ / χ.
        // ∂_i j_c^i = (∂_i χ) h^{ij} j_{cj} + χ (∂_i h^{ij}) j_{cj}
        //             + χ h^{ij} ∂_i j_{cj},
        // and ∂_i h^{ij} = − h^{ik} h^{jl} ∂_i h_{kl}.
        amrex::Real div_j = 0.0;
        FOR (i, j)
        {
            div_j += d1.chi(i) * h_UU[i][j] * j_c[j];
            div_j += chi * h_UU[i][j] * d1.jctrl(j, i);
        }
        FOR (i, j)
        {
            FOR (k, l)
            {
                div_j += -chi * h_UU[i][k] * h_UU[j][l] * d1.h(k, l, i) * j_c[j];
            }
        }
        const amrex::Real inv_chi = 1.0 / chi;
        FOR (i) { div_j += -1.5 * j_cU[i] * d1.chi(i) * inv_chi; }

        amrex::Real j_dot_dalpha = 0.0;
        FOR (i) { j_dot_dalpha += j_cU[i] * d1.lapse(i); }

        // Pump 4-force from the shared RLPumpForce law (same source as the RHS).
        const RLPumpForceDensity force = compute_force(vars, d1, coords, time);

        rhs[c_rho_ctrl] += advec.rho_ctrl() + lapse * K * rho_c -
                           lapse * div_j - 2.0 * j_dot_dalpha - force.f_perp;

        Tensor<1, amrex::Real> f_i;
        f_i[0] = force.f_x;
        f_i[1] = force.f_y;
        f_i[2] = force.f_z;

        const int j_comps[3] = {c_jctrl1, c_jctrl2, c_jctrl3};
        FOR (i)
        {
            amrex::Real shift_term = 0.0;
            // j_{ck} ∂_i β^k, with d1.shift(k, i) = ∂_i β^k.
            FOR (k) { shift_term += j_c[k] * d1.shift(k, i); }
            rhs[j_comps[i]] += advec.jctrl(i) + lapse * K * j_c[i] +
                               shift_term - rho_c * d1.lapse(i) + f_i[i];
        }
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
