#if !defined(COMPLEXSCALARFIELD_HPP_)
#error "This file should only be included through ComplexScalarField.hpp"
#endif

#ifndef COMPLEXSCALARFIELD_IMPL_HPP_
#define COMPLEXSCALARFIELD_IMPL_HPP_

AMREX_GPU_DEVICE AMREX_FORCE_INLINE emtensor_t ComplexScalarField::compute_emtensor(
    const Vars &vars, const D1Vars &d1, const Tensor<2, amrex::Real> &h_UU,
    const Tensor<3, amrex::Real> &chris_ULL) const
{
    (void)chris_ULL;
    const amrex::Real sign = m_sign;
    emtensor_t out;
    out.rho = 0.0;
    out.trS = 0.0;
    FOR (i)
    {
        out.j[i] = 0.0;
    }
    FOR (i, j)
    {
        out.S[i][j] = 0.0;
    }

    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi1  = 0.0;
    amrex::Real dVdphi2  = 0.0;
    m_potential.compute_potential(V_of_phi, dVdphi1, dVdphi2, vars.phi1(),
                                  vars.phi2());

    for (int comp = 0; comp < 2; ++comp)
    {
        const amrex::Real Pi_k = (comp == 0) ? vars.Pi1() : vars.Pi2();

        amrex::Real Vt_k = -Pi_k * Pi_k;
        FOR (i, j)
        {
            const amrex::Real dphi_i =
                (comp == 0) ? d1.phi1(i) : d1.phi2(i);
            const amrex::Real dphi_j =
                (comp == 0) ? d1.phi1(j) : d1.phi2(j);
            Vt_k += vars.chi() * h_UU[i][j] * dphi_i * dphi_j;
        }

        out.rho += Pi_k * Pi_k + 0.5 * Vt_k;
        FOR (i)
        {
            const amrex::Real dphi_i =
                (comp == 0) ? d1.phi1(i) : d1.phi2(i);
            out.j[i] += -Pi_k * dphi_i;
        }
        FOR (i, j)
        {
            const amrex::Real dphi_i =
                (comp == 0) ? d1.phi1(i) : d1.phi2(i);
            const amrex::Real dphi_j =
                (comp == 0) ? d1.phi1(j) : d1.phi2(j);
            out.S[i][j] += -0.5 * vars.h(i, j) * Vt_k / vars.chi() +
                           dphi_i * dphi_j;
        }
    }

    out.rho += V_of_phi;
    FOR (i, j)
    {
        out.S[i][j] -= vars.h(i, j) * V_of_phi / vars.chi();
    }

    // Phantom sign: flip entire T_ab for sign == -1. The field EOM (RHS) is
    // unchanged -- only the gravitational coupling is reversed, giving
    // negative energy density while preserving U(1) charge conservation.
    out.rho *= sign;
    FOR (i)
    {
        out.j[i] *= sign;
    }
    FOR (i, j)
    {
        out.S[i][j] *= sign;
    }
    // S is already sign-flipped, so trace(S) carries the sign; only the
    // explicit -3V term needs the extra factor.
    out.trS =
        vars.chi() * TensorAlgebra::compute_trace(out.S, h_UU) -
        sign * 3.0 * V_of_phi;

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE emtensor_t ComplexScalarField::compute_emtensor(
    const Vars &vars, const D1Vars &d1, const Tensor<2, amrex::Real> &h_UU,
    const Tensor<3, amrex::Real> &chris_ULL, const Coordinates &coords,
    amrex::Real time) const
{
    (void)coords;
    (void)time;
    return compute_emtensor(vars, d1, h_UU, chris_ULL);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void ComplexScalarField::add_matter_rhs(
    const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
    const D1Vars &d1, const D2Vars &d2, const AdvecVars &advec) const
{
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi1  = 0.0;
    amrex::Real dVdphi2  = 0.0;
    m_potential.compute_potential(V_of_phi, dVdphi1, dVdphi2, vars.phi1(),
                                  vars.phi2());

    rhs[c_phi]  = vars.lapse() * vars.Pi1() + advec.phi1();
    rhs[c_Pi]   = vars.lapse() * (vars.K() * vars.Pi1() - dVdphi1) + advec.Pi1();
    rhs[c_phi2] = vars.lapse() * vars.Pi2() + advec.phi2();
    rhs[c_Pi2]  = vars.lapse() * (vars.K() * vars.Pi2() - dVdphi2) + advec.Pi2();

    FOR (i, j)
    {
        rhs[c_Pi] += h_UU[i][j] *
                     (-0.5 * d1.chi(j) * vars.lapse() * d1.phi1(i) +
                      vars.chi() * vars.lapse() * d2.phi1[i][j] +
                      vars.chi() * d1.lapse(i) * d1.phi1(j));
        rhs[c_Pi2] += h_UU[i][j] *
                      (-0.5 * d1.chi(j) * vars.lapse() * d1.phi2(i) +
                       vars.chi() * vars.lapse() * d2.phi2[i][j] +
                       vars.chi() * d1.lapse(i) * d1.phi2(j));
        FOR (k)
        {
            rhs[c_Pi] += -vars.chi() * vars.lapse() * h_UU[i][j] *
                         chris.ULL[k][i][j] * d1.phi1(k);
            rhs[c_Pi2] += -vars.chi() * vars.lapse() * h_UU[i][j] *
                          chris.ULL[k][i][j] * d1.phi2(k);
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void ComplexScalarField::add_matter_rhs(
    const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
    const D1Vars &d1, const D2Vars &d2, const AdvecVars &advec,
    const Coordinates &coords, amrex::Real time) const
{
    add_matter_rhs(rhs, vars, d1, d2, advec);

    if (m_pump.num_sites < 1)
    {
        return;
    }

    const amrex::Real governor = m_pump.governor;

    // ---- Closed-loop PD "trap" controller (k_p > 0) ----------------------
    // Drive the field toward the TARGET soliton at the moving trajectory centre:
    // Phi*(x,t) = amp * gauss(x - centre(t)) * e^{i*arg}, arg = -omega*t + phase,
    // with target momentum Pi* = (1/alpha) d_t Phi*.  Error-proportional, hence
    // self-limiting: confines/transports a coherent soliton along the trajectory
    // rather than dumping energy in (which disperses the lump).
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
            // Trap-target central amplitude: the star's own amplitude
            // (target_amp = bs_phi_c) when set, else the site amplitude.  Aiming
            // at the correct amplitude is essential -- a mismatched (e.g. 10x
            // too small) target makes the trap fight gravity and collapse the star.
            const amrex::Real amp_t =
                (m_pump.target_amp > 0.0) ? m_pump.target_amp : site.amplitude;
            const amrex::Real g   = amp_t * env;
            const amrex::Real arg = -site.frequency * time + site.phase;
            const amrex::Real tphi1 = g * std::cos(arg);
            const amrex::Real tphi2 = g * std::sin(arg);
            const amrex::Real tPi1  = site.frequency * tphi2 * inv_alpha;
            const amrex::Real tPi2  = -site.frequency * tphi1 * inv_alpha;
            const amrex::Real w     = governor * env;
            rhs[c_Pi] +=
                w * (-kp * (vars.phi1() - tphi1) - kd * (vars.Pi1() - tPi1));
            rhs[c_Pi2] +=
                w * (-kp * (vars.phi2() - tphi2) - kd * (vars.Pi2() - tPi2));
        }
        return;
    }

    // ---- Legacy open-loop source pump (k_p <= 0) -------------------------
    // One spotlight per lump.  All sites drive the same complex field
    // components (c_Pi / c_Pi2); they reach different lumps because each
    // envelope is localized at that lump's tracked 3-D centre.  Pi1/Pi2 are
    // driven 90 deg out of phase => local U(1) (Noether-charge) injection.
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
        rhs[c_Pi] += base * std::cos(arg);
        rhs[c_Pi2] += base * std::sin(arg);
    }
}

#endif /* COMPLEXSCALARFIELD_IMPL_HPP_ */
