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
#include "ScalarFieldAdvecVars.hpp"
#include "ScalarFieldD1Vars.hpp"
#include "ScalarFieldD2Vars.hpp"
#include "ScalarFieldVars.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS, total num of components
#include "Tensor.hpp"
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
template <class potential_t = DefaultPotential> class ScalarField
{
  protected:
    potential_t m_potential;
    //! The local copy of the potential

  public:

    //!  Constructor of class ScalarField, inputs are the matter parameters.
    ScalarField() = default;

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

    //! The function which adds in the RHS for the matter field vars,
    //! including the potential
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_matter_rhs(const amrex::CellData<amrex::Real> &rhs, 
                   const Vars &vars, const D1Vars &d1,
                   const D2Vars &d2, const AdvecVars &advec) const;
};

#include "ScalarField.impl.hpp"

#endif /* SCALARFIELD_HPP_ */
