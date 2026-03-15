/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef EXOTICSCALARFIELD_HPP_
#define EXOTICSCALARFIELD_HPP_

#include "CCZ4Geometry.hpp"
#include "Coordinates.hpp"
#include "DefaultPotential.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ScalarFieldAdvecVars.hpp"
#include "ScalarFieldD1Vars.hpp"
#include "ScalarFieldD2Vars.hpp"
#include "ScalarFieldVars.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS, total num of components
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "GRParmParse.hpp"

//!  Calculates the matter type specific elements such as the EMTensor and
//   matter evolution
/*!
     This class is an example of a matter_t object which calculates the
     matter type specific elements for the RHS update and the evaluation
     of the constraints. This includes the Energy Momentum Tensor, and
     the matter evolution terms. In this case, a scalar field,
     the matter elements are phi and (minus) its conjugate momentum, Pi.
     It is templated over a potential function potential_t which the
     user must specify in a class, although a default is provided which
     sets dVdphi and V_of_phi to zero.
     It assumes minimal coupling of the field to gravity.
     \sa MatterCCZ4(), ConstraintsMatter()
*/
template <class potential_t = DefaultPotential> class ExoticScalarField
{
  protected:
    potential_t m_potential;
    double m_support_strength;
    double m_support_ramp_start;
    double m_support_ramp_duration;
    double m_support_causal_speed;
    int m_metric_type;
    double m_throat_radius_A;
    double m_throat_radius_B;
    std::array<double, AMREX_SPACEDIM> m_centerA;
    std::array<double, AMREX_SPACEDIM> m_centerB;
    //! The local copy of the potential

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
    support_schedule(amrex::Real time) const;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    nearest_throat_radius(const Coordinates &coords) const;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
    local_support_strength(const Coordinates &coords, amrex::Real time) const;

  public:

    //!  Constructor of class ExoticScalarField, inputs are the matter parameters.
    ExoticScalarField(potential_t a_potential = potential_t(), double a_support_strength = 1.0)
        : m_potential(a_potential), m_support_strength(a_support_strength),
          m_support_ramp_start(-1.0), m_support_ramp_duration(5.0),
          m_support_causal_speed(0.0), m_metric_type(1),
          m_throat_radius_A(1.0), m_throat_radius_B(1.0), m_centerA{0.0},
          m_centerB{0.0}
    {
        GRParmParse pp;
        if (a_support_strength == 1.0)
            pp.load("wormhole_support_strength", m_support_strength, 1.0);

        pp.load("wormhole_support_ramp_start", m_support_ramp_start, -1.0);
        pp.load("wormhole_support_ramp_duration", m_support_ramp_duration, 5.0);
        pp.load("wormhole_support_causal_speed", m_support_causal_speed, 0.0);
        pp.load("wormhole_metric_type", m_metric_type, 1);
        pp.load("wormhole_throat_radius_A", m_throat_radius_A, 1.0);
        pp.load("wormhole_throat_radius_B", m_throat_radius_B, m_throat_radius_A);

        std::array<double, AMREX_SPACEDIM> default_centerA = {0.0, 0.0, 0.0};
        pp.load("wormhole_centerA", m_centerA, default_centerA);

        std::array<double, AMREX_SPACEDIM> default_centerB = {
            -m_centerA[0], -m_centerA[1], -m_centerA[2]};
        pp.load("wormhole_centerB", m_centerB, default_centerB);
    }

    using Vars      = ScalarFieldVars;
    using D1Vars    = ScalarFieldD1Vars;
    using D2Vars    = ScalarFieldD2Vars;
    using AdvecVars = ScalarFieldAdvecVars;

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, including the potential
    [[nodiscard]]
    AMREX_GPU_DEVICE emtensor_t
    compute_emtensor(const Vars &vars, const D1Vars &d1,
                     const Tensor<2, amrex::Real>
                         &h_UU, //!< the inverse metric (raised indices)
                     const Tensor<3, amrex::Real> &chris_ULL)
        const; //!< the conformal christoffel symbol

    [[nodiscard]]
    AMREX_GPU_DEVICE emtensor_t
    compute_emtensor(const Vars &vars, const D1Vars &d1,
                     const Tensor<2, amrex::Real> &h_UU,
                     const Tensor<3, amrex::Real> &chris_ULL,
                     const Coordinates &coords, amrex::Real time) const;

    //! The function which adds in the RHS for the matter field vars,
    //! including the potential
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_matter_rhs(const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
                   const D1Vars &d1, const D2Vars &d2,
                   const AdvecVars &advec) const;
};

#include "ExoticScalarField.impl.hpp"

#endif /* EXOTICSCALARFIELD_HPP_ */
