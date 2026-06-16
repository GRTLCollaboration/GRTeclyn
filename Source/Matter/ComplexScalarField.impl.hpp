#if !defined(COMPLEXSCALARFIELD_HPP_)
#error "This file should only be included through ComplexScalarField.hpp"
#endif

#ifndef COMPLEXSCALARFIELD_IMPL_HPP_
#define COMPLEXSCALARFIELD_IMPL_HPP_

AMREX_GPU_DEVICE emtensor_t ComplexScalarField::compute_emtensor(
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

AMREX_GPU_DEVICE emtensor_t ComplexScalarField::compute_emtensor(
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

#endif /* COMPLEXSCALARFIELD_IMPL_HPP_ */
