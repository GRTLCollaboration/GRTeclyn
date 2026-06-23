#if !defined(GRTRESNA_INDEPENDENT_SCALARS_HPP_)
#error "This file should only be included through GRTresnaIndependentScalars.hpp"
#endif

#ifndef GRTRESNA_INDEPENDENT_SCALARS_IMPL_HPP_
#define GRTRESNA_INDEPENDENT_SCALARS_IMPL_HPP_

AMREX_GPU_DEVICE AMREX_FORCE_INLINE emtensor_t GRTresnaIndependentScalars::compute_emtensor(
    const Vars &vars, const D1Vars &d1, const Tensor<2, amrex::Real> &h_UU,
    const Tensor<3, amrex::Real> &chris_ULL) const
{
    (void)chris_ULL;
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

    for (int k = 0; k < m_num_fields; ++k)
    {
        const amrex::Real sign = static_cast<amrex::Real>(m_signs[k]);
        amrex::Real Vt_k       = -vars.Pi(k) * vars.Pi(k);
        FOR (i, j)
        {
            Vt_k += vars.chi() * h_UU[i][j] * d1.phi(k, i) * d1.phi(k, j);
        }

        out.rho += sign * (0.5 * vars.Pi(k) * vars.Pi(k) + 0.5 * Vt_k);
        FOR (i)
        {
            out.j[i] += sign * (-d1.phi(k, i) * vars.Pi(k));
        }
        FOR (i, j)
        {
            out.S[i][j] += sign * (-0.5 * vars.h(i, j) * Vt_k / vars.chi() +
                                   d1.phi(k, i) * d1.phi(k, j));
        }
    }

    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi   = 0.0;
    m_potential.compute_potential(V_of_phi, dVdphi, phi_sum(vars));
    out.rho += V_of_phi;
    FOR (i, j)
    {
        out.S[i][j] -= vars.h(i, j) * V_of_phi / vars.chi();
    }
    out.trS =
        vars.chi() * TensorAlgebra::compute_trace(out.S, h_UU) - 3.0 * V_of_phi;

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE emtensor_t GRTresnaIndependentScalars::compute_emtensor(
    const Vars &vars, const D1Vars &d1, const Tensor<2, amrex::Real> &h_UU,
    const Tensor<3, amrex::Real> &chris_ULL, const Coordinates &coords,
    amrex::Real time) const
{
    (void)coords;
    (void)time;
    return compute_emtensor(vars, d1, h_UU, chris_ULL);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void GRTresnaIndependentScalars::add_matter_rhs(
    const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
    const D1Vars &d1, const D2Vars &d2, const AdvecVars &advec) const
{
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi   = 0.0;
    m_potential.compute_potential(V_of_phi, dVdphi, phi_sum(vars));

    for (int k = 0; k < m_num_fields; ++k)
    {
        rhs[c_phi_lump_index(k)] = vars.lapse() * vars.Pi(k) + advec.phi(k);
        rhs[c_Pi_lump_index(k)] =
            vars.lapse() * (vars.K() * vars.Pi(k) - dVdphi) + advec.Pi(k);

        FOR (i, j)
        {
            rhs[c_Pi_lump_index(k)] +=
                h_UU[i][j] *
                (-0.5 * d1.chi(j) * vars.lapse() * d1.phi(k, i) +
                 vars.chi() * vars.lapse() * d2.phi(k)[i][j] +
                 vars.chi() * d1.lapse(i) * d1.phi(k, j));
            FOR (l)
            {
                rhs[c_Pi_lump_index(k)] +=
                    -vars.chi() * vars.lapse() * h_UU[i][j] *
                    chris.ULL[l][i][j] * d1.phi(k, l);
            }
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void GRTresnaIndependentScalars::add_matter_rhs(
    const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
    const D1Vars &d1, const D2Vars &d2, const AdvecVars &advec,
    const Coordinates &coords, amrex::Real time) const
{
    add_matter_rhs(rhs, vars, d1, d2, advec);

    if (m_pump.num_sites < 1)
    {
        return;
    }

    // One spotlight per lump: site s drives lump s's conjugate momentum.
    const amrex::Real governor = m_pump.governor;
    const int n_sites =
        (m_pump.num_sites < m_num_fields) ? m_pump.num_sites : m_num_fields;
    for (int s = 0; s < n_sites; ++s)
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
        rhs[c_Pi_lump_index(s)] += base * std::cos(arg);
    }
}

#endif /* GRTRESNA_INDEPENDENT_SCALARS_IMPL_HPP_ */
