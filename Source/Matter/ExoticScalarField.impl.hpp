/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(EXOTICSCALARFIELD_HPP_)
#error "This file should only be included through ExoticScalarField.hpp"
#endif

#ifndef EXOTICSCALARFIELD_IMPL_HPP_
#define EXOTICSCALARFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t>
AMREX_GPU_DEVICE emtensor_t ExoticScalarField<potential_t>::compute_emtensor(
    const Vars &vars, const D1Vars &d1, const Tensor<2, amrex::Real> &h_UU,
    const Tensor<3, amrex::Real> &chris_ULL) const
{
    const double support_strength = m_support_strength;
    emtensor_t out;

    //    Useful quantity Vt
    amrex::Real Vt = -vars.Pi() * vars.Pi();
    FOR (i, j)
    {
        Vt += vars.chi() * h_UU[i][j] * d1.phi(i) * d1.phi(j);
    }

    // set the potential values
    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi   = 0.0;

    // compute potential and add constributions to EM Tensor
    m_potential.compute_potential(V_of_phi, dVdphi, vars);

    // Calculate components of EM Tensor
    // S = T_ij
    FOR (i, j)
    {
        out.S[i][j] = -support_strength *
                      (-0.5 * vars.h(i, j) * Vt / vars.chi() + d1.phi(i) * d1.phi(j) -
                       vars.h(i, j) * V_of_phi / vars.chi());
    }

    // rho = n^a n^b T_ab
    out.rho = -support_strength * (vars.Pi() * vars.Pi() + 0.5 * Vt + V_of_phi);

    // trS = Tr_S_ij
    // note: out.S is already scaled by support_strength, so compute_trace(out.S)
    // is scaled. normal trS has -3.0 * V, so phantom trS has +3.0 * V.
    out.trS = vars.chi() * TensorAlgebra::compute_trace(out.S, h_UU) +
              support_strength * 3.0 * V_of_phi;

    //    j_i (note lower index) = - n^a T_ai
    FOR (i)
    {
        out.j[i] = -support_strength * (-d1.phi(i) * vars.Pi());
    }

    return out;
}

template <class potential_t>
AMREX_GPU_DEVICE emtensor_t ExoticScalarField<potential_t>::compute_emtensor(
    const Vars &vars, const D1Vars &d1, const Tensor<2, amrex::Real> &h_UU,
    const Tensor<3, amrex::Real> &chris_ULL, const Coordinates &coords,
    amrex::Real time) const
{
    const double support_strength = local_support_strength(coords, time);
    emtensor_t out;

    //    Useful quantity Vt
    amrex::Real Vt = -vars.Pi() * vars.Pi();
    FOR (i, j)
    {
        Vt += vars.chi() * h_UU[i][j] * d1.phi(i) * d1.phi(j);
    }

    // set the potential values
    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi   = 0.0;

    // compute potential and add constributions to EM Tensor
    m_potential.compute_potential(V_of_phi, dVdphi, vars);

    // Calculate components of EM Tensor
    // S = T_ij
    FOR (i, j)
    {
        out.S[i][j] = -support_strength *
                      (-0.5 * vars.h(i, j) * Vt / vars.chi() + d1.phi(i) * d1.phi(j) -
                       vars.h(i, j) * V_of_phi / vars.chi());
    }

    // rho = n^a n^b T_ab
    out.rho = -support_strength * (vars.Pi() * vars.Pi() + 0.5 * Vt + V_of_phi);

    // trS = Tr_S_ij
    // note: out.S is already scaled by support_strength, so compute_trace(out.S)
    // is scaled.
    // normal trS has -3.0 * V, so phantom trS has +3.0 * V.
    out.trS = vars.chi() * TensorAlgebra::compute_trace(out.S, h_UU) +
              support_strength * 3.0 * V_of_phi;

    //    j_i (note lower index) = - n^a T_ai
    FOR (i)
    {
        out.j[i] = -support_strength * (-d1.phi(i) * vars.Pi());
    }

    return out;
}

template <class potential_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
ExoticScalarField<potential_t>::support_schedule(amrex::Real time) const
{
    if (m_support_ramp_start < 0.0 || time <= m_support_ramp_start)
        return 1.0;

    if (m_support_ramp_duration <= 0.0 ||
        time >= m_support_ramp_start + m_support_ramp_duration)
        return 0.0;

    const double x = (time - m_support_ramp_start) / m_support_ramp_duration;
    return 0.5 * (1.0 + cos(M_PI * x));
}

template <class potential_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
ExoticScalarField<potential_t>::nearest_throat_radius(
    const Coordinates &coords) const
{
    const amrex::Real dxA = coords.x - m_centerA[0];
    const amrex::Real dyA = coords.y - m_centerA[1];
    const amrex::Real dzA = coords.z - m_centerA[2];
    const amrex::Real rA =
        std::sqrt(dxA * dxA + dyA * dyA + dzA * dzA);
    const amrex::Real throatA = amrex::Math::abs(
        rA - amrex::Real(0.5 * m_throat_radius_A));

    if (m_metric_type == 1)
        return throatA;

    const amrex::Real dxB = coords.x - m_centerB[0];
    const amrex::Real dyB = coords.y - m_centerB[1];
    const amrex::Real dzB = coords.z - m_centerB[2];
    const amrex::Real rB =
        std::sqrt(dxB * dxB + dyB * dyB + dzB * dzB);
    const amrex::Real throatB = amrex::Math::abs(
        rB - amrex::Real(0.5 * m_throat_radius_B));
    return amrex::min(throatA, throatB);
}

template <class potential_t>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
ExoticScalarField<potential_t>::local_support_strength(
    const Coordinates &coords, amrex::Real time) const
{
    if (m_support_causal_speed <= 0.0)
        return m_support_strength * support_schedule(time);

    const amrex::Real r = nearest_throat_radius(coords);
    const amrex::Real retarded_time = time - r / m_support_causal_speed;
    return m_support_strength * support_schedule(retarded_time);
}

// Adds in the RHS for the matter vars
template <class potential_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ExoticScalarField<potential_t>::add_matter_rhs(
    const amrex::CellData<amrex::Real> &rhs, const Vars &vars, const D1Vars &d1,
    const D2Vars &d2, const AdvecVars &advec) const
{
    // call the function for the rhs excluding the potential
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    // set the potential values
    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi   = 0.0;
    m_potential.compute_potential(V_of_phi, dVdphi, vars);

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs[c_phi] = vars.lapse() * vars.Pi() + advec.phi();

    rhs[c_Pi] = vars.lapse() * (vars.K() * vars.Pi() - dVdphi) + advec.Pi();

    FOR (i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs[c_Pi] += h_UU[i][j] * (-0.5 * d1.chi(j) * vars.lapse() * d1.phi(i) +
                                   vars.chi() * vars.lapse() * d2.phi[i][j] +
                                   vars.chi() * d1.lapse(i) * d1.phi(j));
        FOR (k)
        {
            rhs[c_Pi] += -vars.chi() * vars.lapse() * h_UU[i][j] *
                         chris.ULL[k][i][j] * d1.phi(k);
        }
    }
}

#endif /* EXOTICSCALARFIELD_IMPL_HPP_ */
