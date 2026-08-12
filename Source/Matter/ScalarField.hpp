/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELD_HPP_
#define SCALARFIELD_HPP_

#include "CCZ4Geometry.hpp"
#include "DefaultPotential.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ScalarFieldVars.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS, total num of components
#include "TensorAlgebra.hpp"

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
template <class potential_t = DefaultPotential,
          class deriv_t     = FourthOrderDerivatives>
class ScalarField
{
  protected:
    potential_t m_potential;
    //! The local copy of the potential

  public:

    //!  Constructor of class ScalarField, inputs are the matter parameters.
    ScalarField() = default;

    AMREX_GPU_HOST_DEVICE
        AMREX_FORCE_INLINE explicit ScalarField(potential_t a_potential)
        : m_potential(a_potential)
    {
    }

    using Vars = ScalarFieldVars;

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, including the potential
    [[nodiscard]]
    AMREX_GPU_DEVICE emtensor_t compute_emtensor(
        const int ix, const int iy, const int iz, //!< grid indicies
        const amrex::Array4<const amrex::Real>
            &state,             //!< the current value of state variables
        const deriv_t &a_deriv, //!< the object that calculates the derivative
        const Tensor::Rank2 &h_UU) //!< the inverse metric (raised indices)
        const;

    // ! The function which adds in the RHS for the matter field vars,
    // ! including the potential
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void add_matter_rhs(
        const int ix, const int iy, const int iz, //!< grid indicies
        const amrex::Array4<amrex::Real>
            &rhs_state, //!< the next value of state variables (rhs update)
        const amrex::Array4<const amrex::Real>
            &state, //!< the current value of state variables
        const deriv_t &a_deriv)
        const; //!< the object for calculating derivatives
};

#include "ScalarField.impl.hpp"

#endif /* SCALARFIELD_HPP_ */
