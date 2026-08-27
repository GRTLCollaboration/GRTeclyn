/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CONSTRAINTSWITHMATTER_HPP_
#define CONSTRAINTSWITHMATTER_HPP_

#include "CCZ4Geometry.hpp"
#include "Constraints.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include <array>

//!  Calculates the Hamiltonian and Momentum constraints with matter fields
/*!
     The class calculates the Hamiltonian and Momentum constraints at each point
   in a box. It inherits from the Constraints class which calculates the
   constraints without the matter terms. It adds in the matter terms for a given
   matter class matter_t, which must provide the stress-energy sources with
   their gravitational coupling already applied. For an example of a matter_t
   class see ScalarField. \sa Constraints(), ScalarField()
*/
template <class matter_t> class ConstraintsWithMatter : public Constraints
{
  public:

    //! Constructor of class ConstraintsWithMatter
    /*!
        Can specify the vars of the constraint vars instead of using the
        hardcoded ones.
    */
    ConstraintsWithMatter(amrex::Real dx, int a_c_Ham, const Interval &a_c_Moms,
                          int a_c_Ham_abs_terms              = -1,
                          const Interval &a_c_Moms_abs_terms = Interval());

    //! The compute member which calculates the constraints at each point in the
    //! box
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &constraints,
               const amrex::Array4<amrex::Real const> &state) const;

    static void set_up(int a_state_index, bool a_calc_mom_norm = false);

    // Has signature of DeriveFuncMF so that it can be stored in the derive_lst
    static void compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
                           const amrex::MultiFab &src_mf,
                           const amrex::Geometry &geomdata,
                           amrex::Real /*time*/, const int * /*bcrec*/,
                           int /*level*/);

  protected:
    matter_t m_matter; //!< The matter object, e.g. a scalar field
};

#include "ConstraintsWithMatter.impl.hpp"

#endif /* CONSTRAINTSWITHMATTER_HPP_ */
