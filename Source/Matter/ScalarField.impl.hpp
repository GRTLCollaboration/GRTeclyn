/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(SCALARFIELD_HPP_)
#error "This file should only be included through ScalarField.hpp"
#endif

#ifndef SCALARFIELD_IMPL_HPP_
#define SCALARFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t, class deriv_t>
AMREX_GPU_DEVICE emtensor_t ScalarField<potential_t, deriv_t>::compute_emtensor(
    const Vars &vars, const D1Vars &d1_scalar, const TensorArray::Rank2 &h_UU,
    const TensorArray::Rank3 &chris_ULL) const
{
    emtensor_t out;

    //    Useful quantity Vt
    amrex::Real Vt = -vars.Pi() * vars.Pi();
    FOR (i, j)
    {
        Vt += vars.chi() * h_UU(i, j) * d1_scalar.phi(i) * d1_scalar.phi(j);
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
        out.S(i, j) = -0.5 * vars.h(i, j) * Vt / vars.chi() +
                      d1_scalar.phi(i) * d1_scalar.phi(j) -
                      vars.h(i, j) * V_of_phi / vars.chi();
    }

    // rho = n^a n^b T_ab
    out.rho = vars.Pi() * vars.Pi() + 0.5 * Vt + V_of_phi;

    // trS = Tr_S_ij
    out.trS =
        vars.chi() * TensorAlgebra::compute_trace(out.S, h_UU) - 3.0 * V_of_phi;

    //    j_i (note lower index) = - n^a T_ai
    FOR (i)
    {
        out.j(i) = -d1_scalar.phi(i) * vars.Pi();
    }

    return out;
}

// Adds in the RHS for the matter vars
template <class potential_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ScalarField<potential_t, deriv_t>::add_matter_rhs(
    const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
    const TensorArray::Rank1 &d1_chi, const TensorArray::Rank1 &d1_lapse,
    const amrex::Array2D<amrex::Real, 0, UNIQUE_IDX - 1, 0, AMREX_SPACEDIM - 1>
        &d1_h,
    const D1Vars &d1_scalar, const D2Vars &d2_scalar,
    const AdvecVars &advec) const
{
    // call the function for the rhs excluding the potential
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

    // set the potential values
    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi   = 0.0;
    m_potential.compute_potential(V_of_phi, dVdphi, vars);

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs[c_phi] = vars.lapse() * vars.Pi() + advec.phi;

    rhs[c_Pi] = vars.lapse() * (vars.K() * vars.Pi() - dVdphi) + advec.Pi;

    FOR (i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs[c_Pi] +=
            h_UU(i, j) *
            (-0.5 * d1_chi(j) * vars.lapse() * d1_scalar.phi(i) +
             vars.chi() * vars.lapse() * d2_scalar.phi(VAR_IDX0(i, j)) +
             vars.chi() * d1_lapse(i) * d1_scalar.phi(j));
        FOR (k)
        {
            rhs[c_Pi] += -vars.chi() * vars.lapse() * h_UU(i, j) *
                         chris.ULL(k, i, j) * d1_scalar.phi(k);
        }
    }
}

#endif /* SCALARFIELD_IMPL_HPP_ */
