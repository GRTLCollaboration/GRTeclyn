/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef EMTENSOR_HPP
#define EMTENSOR_HPP

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.hpp"
#include "MatterCCZ4.hpp"
#include "simd.hpp"

//! Calculates the EM tensor and then saves the ones specified in the
//! constructor on the grid
template <class matter_t> class EMTensor
{
  public:
    template <class data_t>
    using Vars = typename MatterCCZ4<matter_t>::template Vars<data_t>;

    //! Constructor
    EMTensor(const matter_t &a_matter, const double dx, const int a_c_rho = -1,
             const Interval a_c_Si  = Interval(),
             const Interval a_c_Sij = Interval());

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute(int i, int j, int k, const amrex::Array4<data_t> &out_mf,
            const amrex::Array4<const data_t> &in_mf) const;

  protected:
    const matter_t &m_matter;
    FourthOrderDerivatives m_deriv;
    const int m_c_rho;      // var enum for the energy density
    const Interval m_c_Si;  // Interval of var enums for the momentum density
    const Interval m_c_Sij; // Interval of var enums for the spatial
                            // stress-energy density
};

#include "EMTensor.impl.hpp"

#endif /* EMTENSOR_HPP */
