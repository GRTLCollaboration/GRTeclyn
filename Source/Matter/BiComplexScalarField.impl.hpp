#if !defined(BICOMPLEXSCALARFIELD_HPP_)
#error "This file should only be included through BiComplexScalarField.hpp"
#endif

#ifndef BICOMPLEXSCALARFIELD_IMPL_HPP_
#define BICOMPLEXSCALARFIELD_IMPL_HPP_

namespace
{
//! Accumulate one complex field's (signed) stress-energy into out.  Mirrors
//! the single-field ComplexScalarField math; ``fs`` is the field's
//! gravitational sign (+1 canonical, -1 phantom).  Returns fs * V (the signed
//! potential, needed for the trace correction).
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real accumulate_complex_field(
    emtensor_t &out, amrex::Real fs, amrex::Real Pi1, amrex::Real Pi2,
    const Tensor<1, amrex::Real> &d1phi1, const Tensor<1, amrex::Real> &d1phi2,
    amrex::Real V_of_phi, amrex::Real chi, const BiComplexScalarFieldVars &vars,
    const Tensor<2, amrex::Real> &h_UU)
{
    for (int comp = 0; comp < 2; ++comp)
    {
        const amrex::Real Pi_k = (comp == 0) ? Pi1 : Pi2;

        amrex::Real Vt_k = -Pi_k * Pi_k;
        FOR (i, j)
        {
            const amrex::Real dphi_i = (comp == 0) ? d1phi1[i] : d1phi2[i];
            const amrex::Real dphi_j = (comp == 0) ? d1phi1[j] : d1phi2[j];
            Vt_k += chi * h_UU[i][j] * dphi_i * dphi_j;
        }

        out.rho += fs * (Pi_k * Pi_k + 0.5 * Vt_k);
        FOR (i)
        {
            const amrex::Real dphi_i = (comp == 0) ? d1phi1[i] : d1phi2[i];
            out.j[i] += fs * (-Pi_k * dphi_i);
        }
        FOR (i, j)
        {
            const amrex::Real dphi_i = (comp == 0) ? d1phi1[i] : d1phi2[i];
            const amrex::Real dphi_j = (comp == 0) ? d1phi1[j] : d1phi2[j];
            out.S[i][j] += fs * (-0.5 * vars.h(i, j) * Vt_k / chi +
                                 dphi_i * dphi_j);
        }
    }

    out.rho += fs * V_of_phi;
    FOR (i, j)
    {
        out.S[i][j] -= fs * vars.h(i, j) * V_of_phi / chi;
    }
    return fs * V_of_phi;
}
} // namespace

AMREX_GPU_DEVICE AMREX_FORCE_INLINE emtensor_t
BiComplexScalarField::compute_emtensor(
    const Vars &vars, const D1Vars &d1, const Tensor<2, amrex::Real> &h_UU,
    const Tensor<3, amrex::Real> &chris_ULL) const
{
    (void)chris_ULL;
    emtensor_t out;
    out.rho = 0.0;
    out.trS = 0.0;
    FOR (i) { out.j[i] = 0.0; }
    FOR (i, j) { out.S[i][j] = 0.0; }

    // Canonical (Phi+) field, gravitational sign +1.
    amrex::Real Vp = 0.0, dVp1 = 0.0, dVp2 = 0.0;
    m_potential.compute_potential(Vp, dVp1, dVp2, vars.phi1p(), vars.phi2p());
    amrex::Real sumV = accumulate_complex_field(
        out, +1.0, vars.Pi1p(), vars.Pi2p(), d1.m_phi1p_d1, d1.m_phi2p_d1, Vp,
        vars.chi(), vars, h_UU);

    // Phantom (Phi-) field, gravitational sign -1.
    amrex::Real Vm = 0.0, dVm1 = 0.0, dVm2 = 0.0;
    m_potential.compute_potential(Vm, dVm1, dVm2, vars.phi1m(), vars.phi2m());
    sumV += accumulate_complex_field(out, -1.0, vars.Pi1m(), vars.Pi2m(),
                                     d1.m_phi1m_d1, d1.m_phi2m_d1, Vm,
                                     vars.chi(), vars, h_UU);

    out.trS = vars.chi() * TensorAlgebra::compute_trace(out.S, h_UU) -
              3.0 * sumV;
    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE emtensor_t
BiComplexScalarField::compute_emtensor(
    const Vars &vars, const D1Vars &d1, const Tensor<2, amrex::Real> &h_UU,
    const Tensor<3, amrex::Real> &chris_ULL, const Coordinates &coords,
    amrex::Real time) const
{
    (void)coords;
    (void)time;
    return compute_emtensor(vars, d1, h_UU, chris_ULL);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void BiComplexScalarField::add_matter_rhs(
    const amrex::CellData<amrex::Real> &rhs, const Vars &vars, const D1Vars &d1,
    const D2Vars &d2, const AdvecVars &advec) const
{
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    // Potential derivatives for each field (the field EOM is sign-independent;
    // only the gravitational coupling is reversed, so both fields obey the same
    // Klein-Gordon RHS while preserving their own U(1) charge).
    amrex::Real Vp = 0.0, dVp1 = 0.0, dVp2 = 0.0;
    m_potential.compute_potential(Vp, dVp1, dVp2, vars.phi1p(), vars.phi2p());
    amrex::Real Vm = 0.0, dVm1 = 0.0, dVm2 = 0.0;
    m_potential.compute_potential(Vm, dVm1, dVm2, vars.phi1m(), vars.phi2m());

    // Canonical field Phi+.
    rhs[c_phi]  = vars.lapse() * vars.Pi1p() + advec.phi1p();
    rhs[c_Pi]   = vars.lapse() * (vars.K() * vars.Pi1p() - dVp1) + advec.Pi1p();
    rhs[c_phi2] = vars.lapse() * vars.Pi2p() + advec.phi2p();
    rhs[c_Pi2]  = vars.lapse() * (vars.K() * vars.Pi2p() - dVp2) + advec.Pi2p();

    // Phantom field Phi-.
    rhs[c_phi_m]  = vars.lapse() * vars.Pi1m() + advec.phi1m();
    rhs[c_Pi_m]   = vars.lapse() * (vars.K() * vars.Pi1m() - dVm1) + advec.Pi1m();
    rhs[c_phi2_m] = vars.lapse() * vars.Pi2m() + advec.phi2m();
    rhs[c_Pi2_m]  = vars.lapse() * (vars.K() * vars.Pi2m() - dVm2) + advec.Pi2m();

    FOR (i, j)
    {
        rhs[c_Pi] += h_UU[i][j] *
                     (-0.5 * d1.chi(j) * vars.lapse() * d1.phi1p(i) +
                      vars.chi() * vars.lapse() * d2.phi1p[i][j] +
                      vars.chi() * d1.lapse(i) * d1.phi1p(j));
        rhs[c_Pi2] += h_UU[i][j] *
                      (-0.5 * d1.chi(j) * vars.lapse() * d1.phi2p(i) +
                       vars.chi() * vars.lapse() * d2.phi2p[i][j] +
                       vars.chi() * d1.lapse(i) * d1.phi2p(j));
        rhs[c_Pi_m] += h_UU[i][j] *
                       (-0.5 * d1.chi(j) * vars.lapse() * d1.phi1m(i) +
                        vars.chi() * vars.lapse() * d2.phi1m[i][j] +
                        vars.chi() * d1.lapse(i) * d1.phi1m(j));
        rhs[c_Pi2_m] += h_UU[i][j] *
                        (-0.5 * d1.chi(j) * vars.lapse() * d1.phi2m(i) +
                         vars.chi() * vars.lapse() * d2.phi2m[i][j] +
                         vars.chi() * d1.lapse(i) * d1.phi2m(j));
        FOR (k)
        {
            const amrex::Real cc =
                -vars.chi() * vars.lapse() * h_UU[i][j] * chris.ULL[k][i][j];
            rhs[c_Pi]    += cc * d1.phi1p(k);
            rhs[c_Pi2]   += cc * d1.phi2p(k);
            rhs[c_Pi_m]  += cc * d1.phi1m(k);
            rhs[c_Pi2_m] += cc * d1.phi2m(k);
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void BiComplexScalarField::add_matter_rhs(
    const amrex::CellData<amrex::Real> &rhs, const Vars &vars, const D1Vars &d1,
    const D2Vars &d2, const AdvecVars &advec, const Coordinates &coords,
    amrex::Real time) const
{
    add_matter_rhs(rhs, vars, d1, d2, advec);

    if (m_pump.num_sites < 1)
    {
        return;
    }

    const amrex::Real governor = m_pump.governor;

    // ---- Closed-loop PD "trap" controller (k_p > 0) ----------------------
    // Each spotlight drives its field toward the TARGET soliton at the moving
    // trajectory centre:  Phi*(x,t) = amp * gauss(x - centre(t)) * e^{i*arg},
    // arg = -omega*t + phase, with target momentum Pi* = (1/alpha) d_t Phi*.
    // The drive is proportional to the ERROR (field - target), so it is
    // self-limiting: it adds matter where the lump is deficient AND removes it
    // where it strayed, transporting a coherent soliton along the trajectory
    // instead of bulldozing energy in.  Routed to the canonical (field_sign>=0)
    // or phantom field by sign.
    if (m_pump.k_p > 0.0)
    {
        const amrex::Real kp        = m_pump.k_p;
        const amrex::Real kd        = m_pump.k_d;
        const amrex::Real inv_alpha = 1.0 / vars.lapse();
        const amrex::Real tw =
            (m_pump.target_width > 0.0) ? m_pump.target_width : m_pump.width;
        for (int s = 0; s < m_pump.num_sites; ++s)
        {
            const auto &site = m_pump.sites[s];
            if (site.amplitude <= 0.0)
            {
                continue;
            }
            const amrex::Real env = RLRuntime::compute_site_envelope(
                coords.x, coords.y, coords.z, site, tw, m_pump.target_profile);
            if (env < 1.0e-8)
            {
                continue;
            }
            const amrex::Real g   = site.amplitude * env; // target |phi| here
            const amrex::Real arg = -site.frequency * time + site.phase;
            const amrex::Real cc  = std::cos(arg);
            const amrex::Real ss  = std::sin(arg);
            const amrex::Real tphi1 = g * cc;
            const amrex::Real tphi2 = g * ss;
            // Pi* = (1/alpha) d_t Phi* for Phi* = g e^{i arg}, d(arg)/dt=-omega.
            const amrex::Real tPi1 = site.frequency * tphi2 * inv_alpha;
            const amrex::Real tPi2 = -site.frequency * tphi1 * inv_alpha;
            const amrex::Real w = governor * env; // localization window
            if (site.field_sign >= 0)
            {
                rhs[c_Pi] += w * (-kp * (vars.phi1p() - tphi1) -
                                  kd * (vars.Pi1p() - tPi1));
                rhs[c_Pi2] += w * (-kp * (vars.phi2p() - tphi2) -
                                   kd * (vars.Pi2p() - tPi2));
            }
            else
            {
                rhs[c_Pi_m] += w * (-kp * (vars.phi1m() - tphi1) -
                                    kd * (vars.Pi1m() - tPi1));
                rhs[c_Pi2_m] += w * (-kp * (vars.phi2m() - tphi2) -
                                     kd * (vars.Pi2m() - tPi2));
            }
        }
        return;
    }

    // ---- Legacy open-loop source pump (k_p <= 0) -------------------------
    // One spotlight per lump.  A site drives the CANONICAL field (c_Pi/c_Pi2)
    // when field_sign >= 0, otherwise the PHANTOM field (c_Pi_m/c_Pi2_m).
    // Pi1/Pi2 are driven 90 deg out of phase => local U(1) charge injection.
    for (int s = 0; s < m_pump.num_sites; ++s)
    {
        const amrex::Real base = RLRuntime::compute_site_base(
            coords.x, coords.y, coords.z, m_pump.sites[s], m_pump.width,
            governor);
        if (base <= 0.0)
        {
            continue;
        }
        const amrex::Real arg =
            m_pump.sites[s].frequency * time + m_pump.sites[s].phase;
        const amrex::Real dc = base * std::cos(arg);
        const amrex::Real ds = base * std::sin(arg);
        if (m_pump.sites[s].field_sign >= 0)
        {
            rhs[c_Pi]  += dc;
            rhs[c_Pi2] += ds;
        }
        else
        {
            rhs[c_Pi_m]  += dc;
            rhs[c_Pi2_m] += ds;
        }
    }
}

#endif /* BICOMPLEXSCALARFIELD_IMPL_HPP_ */
