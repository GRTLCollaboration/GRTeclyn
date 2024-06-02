/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MATTERWEYL4_HPP_
#define MATTERWEYL4_HPP_

#include "MatterCCZ4.hpp"
#include "Weyl4.hpp"

//!  Calculates the Weyl4 scalar for spacetimes with matter content
/*!
   This class calculates the Weyl4 scalar real and im parts. It inherits from
   the Weyl4 class and adds in the matter terms as appropriate depending on the
   formulation
*/
template <class matter_t> class MatterWeyl4 : public Weyl4
{
  public:
    template <class data_t>
    using Vars = typename MatterCCZ4<matter_t>::template Vars<data_t>;

    //! Constructor
    MatterWeyl4(matter_t a_matter,
                const std::array<double, AMREX_SPACEDIM> a_center,
                const double a_dx, const int a_dcomp,
                const int a_formulation = CCZ4RHS<>::USE_CCZ4,
                double a_G_Newton       = 1.0)
        : Weyl4(a_center, a_dx, a_formulation), m_matter(a_matter),
          m_dcomp(a_dcomp), m_G_Newton(a_G_Newton)
    {
    }

    //! The compute member which calculates the wave quantities at each point on
    //! the grid
    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute(int i, int j, int k, const amrex::Array4<data_t> &rhs,
            const amrex::Array4<data_t const> &state) const;

  protected:
    matter_t m_matter;       //!< The matter object, e.g. a scalar field
    const int m_dcomp;       //!< index for storing the results of compute
    const double m_G_Newton; //!< Newton's constant, set to one by default

    //! Add matter terms to electric and magnetic parts
    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_matter_EB(EBFields_t<data_t> &eb_fields, const Vars<data_t> &vars,
                  const Vars<Tensor<1, data_t>> &d1,
                  const Tensor<3, data_t> &epsilon3_LUU,
                  const Tensor<2, data_t> &h_UU,
                  const chris_t<data_t> &chris) const;
};

#include "MatterWeyl4.impl.hpp"

#endif /* MATTERWEYL4_HPP_ */
