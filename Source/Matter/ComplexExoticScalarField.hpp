/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef COMPLEXEXOTICSCALARFIELD_HPP_
#define COMPLEXEXOTICSCALARFIELD_HPP_

#include "CCZ4Geometry.hpp"
#include "ComplexScalarFieldAdvecVars.hpp"
#include "ComplexScalarFieldD1Vars.hpp"
#include "ComplexScalarFieldD2Vars.hpp"
#include "ComplexScalarFieldVars.hpp"
#include "Coordinates.hpp"
#include "DefaultPotential.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRParmParse.hpp"
#include "RLMatterPumpParams.hpp"
#include "RLPumpForce.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

//! Exotic (phantom) complex scalar field stored as two real components.
/*!
    The complex field Phi = phi1 + i phi2 has an axisymmetric modulus
    |Phi|^2 = phi1^2 + phi2^2 and carries angular momentum through its phase
    winding, so it can describe a smoothly rotating wormhole without the
    azimuthal density lobes (and outward dispersal) of a single real scalar.

    Energetically each real component contributes like an ExoticScalarField
    component; the overall ``-support_strength`` factor implements the phantom
    (null-energy-condition-violating) sign needed to hold the throat open. The
    potential is templated and evaluated on the full complex modulus via the
    coupled ``compute_potential(V, dV1, dV2, phi1, phi2)`` interface so that a
    self-interacting Q-ball potential V(|Phi|^2) is treated exactly and the U(1)
    (Noether-charge) symmetry is preserved.
*/
template <class potential_t = DefaultPotential> class ComplexExoticScalarField
{
  protected:
    potential_t m_potential;
    double m_support_strength;
    RLMatterPumpParams m_pump{};

  public:
    ComplexExoticScalarField(potential_t a_potential = potential_t(),
                             double a_support_strength = 1.0,
                             RLMatterPumpParams a_pump = {})
        : m_potential(a_potential), m_support_strength(a_support_strength),
          m_pump(a_pump)
    {
        GRParmParse pp;
        if (a_support_strength == 1.0)
            pp.load("wormhole_support_strength", m_support_strength, 1.0);
    }

    using Vars      = ComplexScalarFieldVars;
    using D1Vars    = ComplexScalarFieldD1Vars;
    using D2Vars    = ComplexScalarFieldD2Vars;
    using AdvecVars = ComplexScalarFieldAdvecVars;

    [[nodiscard]]
    AMREX_GPU_DEVICE emtensor_t
    compute_emtensor(const Vars &vars, const D1Vars &d1,
                     const Tensor<2, amrex::Real> &h_UU,
                     const Tensor<3, amrex::Real> &chris_ULL) const;

    [[nodiscard]]
    AMREX_GPU_DEVICE emtensor_t
    compute_emtensor(const Vars &vars, const D1Vars &d1,
                     const Tensor<2, amrex::Real> &h_UU,
                     const Tensor<3, amrex::Real> &chris_ULL,
                     const Coordinates &coords, amrex::Real time) const;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_matter_rhs(const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
                   const D1Vars &d1, const D2Vars &d2,
                   const AdvecVars &advec) const;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_matter_rhs(const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
                   const D1Vars &d1, const D2Vars &d2,
                   const AdvecVars &advec, const Coordinates &coords,
                   amrex::Real time) const;
};

#include "ComplexExoticScalarField.impl.hpp"

#endif /* COMPLEXEXOTICSCALARFIELD_HPP_ */
