/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef WEYL4WITHMATTER_HPP_
#define WEYL4WITHMATTER_HPP_

#include "CCZ4RHSWithMatter.hpp"
#include "Weyl4.hpp"

//!  Calculates the Weyl4 scalar for spacetimes with matter content
/*!
   This class calculates the Weyl4 scalar real and im parts. It inherits from
   the Weyl4 class and adds in the matter terms as appropriate depending on the
   formulation
*/
template <class matter_t> class Weyl4WithMatter : public Weyl4
{
  public:
    //! Constructor
    Weyl4WithMatter(const std::array<double, AMREX_SPACEDIM> a_center,
                    const double a_dx, const int a_dcomp,
                    const int a_formulation = CCZ4RHS<>::USE_CCZ4,
                    double a_G_Newton       = 1.0)
        : Weyl4(a_center, a_dx, a_dcomp, a_formulation), m_dcomp(a_dcomp),
          m_G_Newton(a_G_Newton)
    {
    }

    //! The compute member which calculates the wave quantities at each point on
    //! the grid
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &weyl_scalars,
               const amrex::Array4<amrex::Real const> &state) const;

    static void set_up(int a_state_index);

    // Has signature of DeriveFuncMF so that it can be stored in the derive_lst
    static void compute_mf(amrex::MultiFab &out_mf, int out_comp, int ncomp,
                           const amrex::MultiFab &src_mf,
                           const amrex::Geometry &geomdata,
                           amrex::Real /*time*/, const int * /*bcrec*/,
                           int /*level*/);

  protected:

    matter_t m_matter;
    int m_dcomp;       //!< index for storing the results of compute
    double m_G_Newton; //!< Newton's constant, set to one by default

    //! Add matter terms to electric and magnetic parts
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_matter_EB(EBFields_t &eb_fields, const typename matter_t::Vars &vars,
                  const typename matter_t::D1Vars &d1,
                  const amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3> &epsilon3_LUU,
                  const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &h_UU,
                  const chris_array_t &chris) const;
};

#include "Weyl4WithMatter.impl.hpp"

#endif /* WEYL4WITHMATTER_HPP_ */
