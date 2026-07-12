/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(COMPLEXEXOTICSCALARFIELD_HPP_)
#error "This file should only be included through ComplexExoticScalarField.hpp"
#endif

#ifndef COMPLEXEXOTICSCALARFIELD_IMPL_HPP_
#define COMPLEXEXOTICSCALARFIELD_IMPL_HPP_

template <class potential_t>
AMREX_GPU_DEVICE emtensor_t
ComplexExoticScalarField<potential_t>::compute_emtensor(
    const Vars &vars, const D1Vars &d1, const Tensor<2, amrex::Real> &h_UU,
    const Tensor<3, amrex::Real> &chris_ULL) const
{
    const double support_strength = m_support_strength;
    emtensor_t out;

    // Useful quantity Vt per component
    amrex::Real Vt1 = -vars.Pi1() * vars.Pi1();
    amrex::Real Vt2 = -vars.Pi2() * vars.Pi2();
    FOR (i, j)
    {
        Vt1 += vars.chi() * h_UU[i][j] * d1.phi1(i) * d1.phi1(j);
        Vt2 += vars.chi() * h_UU[i][j] * d1.phi2(i) * d1.phi2(j);
    }
    const amrex::Real Vt = Vt1 + Vt2;

    // Potential, evaluated on the full complex modulus (coupled).  This is
    // exact for a self-interacting potential V(|Phi|^2) = V(phi1^2 + phi2^2);
    // a per-component evaluation would break the U(1) symmetry (and destroy the
    // conserved Noether charge) for the quartic/sextic Q-ball terms.
    amrex::Real V_of_phi = 0.0, dV1 = 0.0, dV2 = 0.0;
    m_potential.compute_potential(V_of_phi, dV1, dV2, vars.phi1(), vars.phi2());

    // S = T_ij (sum of both components, phantom-flipped)
    FOR (i, j)
    {
        out.S[i][j] =
            -support_strength *
            (-0.5 * vars.h(i, j) * Vt / vars.chi() +
             d1.phi1(i) * d1.phi1(j) + d1.phi2(i) * d1.phi2(j) -
             vars.h(i, j) * V_of_phi / vars.chi());
    }

    // rho = n^a n^b T_ab
    out.rho = -support_strength * (vars.Pi1() * vars.Pi1() +
                                   vars.Pi2() * vars.Pi2() + 0.5 * Vt +
                                   V_of_phi);

    // trS: out.S is already scaled by support_strength; restore the +3V piece.
    out.trS = vars.chi() * TensorAlgebra::compute_trace(out.S, h_UU) +
              support_strength * 3.0 * V_of_phi;

    // j_i (lower index) = - n^a T_ai
    FOR (i)
    {
        out.j[i] = -support_strength *
                   (-(d1.phi1(i) * vars.Pi1() + d1.phi2(i) * vars.Pi2()));
    }

    return out;
}

template <class potential_t>
AMREX_GPU_DEVICE emtensor_t
ComplexExoticScalarField<potential_t>::compute_emtensor(
    const Vars &vars, const D1Vars &d1, const Tensor<2, amrex::Real> &h_UU,
    const Tensor<3, amrex::Real> &chris_ULL, const Coordinates & /*coords*/,
    amrex::Real /*time*/) const
{
    return compute_emtensor(vars, d1, h_UU, chris_ULL);
}

template <class potential_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ComplexExoticScalarField<potential_t>::add_matter_rhs(
    const amrex::CellData<amrex::Real> &rhs, const Vars &vars, const D1Vars &d1,
    const D2Vars &d2, const AdvecVars &advec) const
{
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    amrex::Real V_of_phi = 0.0, dV1 = 0.0, dV2 = 0.0;
    m_potential.compute_potential(V_of_phi, dV1, dV2, vars.phi1(), vars.phi2());

    // Component 1: Phi real part
    rhs[c_phi] = vars.lapse() * vars.Pi1() + advec.phi1();
    rhs[c_Pi]  = vars.lapse() * (vars.K() * vars.Pi1() - dV1) + advec.Pi1();

    // Component 2: Phi imaginary part
    rhs[c_phi2] = vars.lapse() * vars.Pi2() + advec.phi2();
    rhs[c_Pi2]  = vars.lapse() * (vars.K() * vars.Pi2() - dV2) + advec.Pi2();

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

#endif /* COMPLEXEXOTICSCALARFIELD_IMPL_HPP_ */
