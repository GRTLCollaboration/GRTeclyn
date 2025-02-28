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

// AMReX Includes
#include <AMReX_MultiFab.H>

//! Calculates the EM tensor and then saves the ones specified in the
//! constructor on the grid
template <class matter_t> class EMTensor
{
  public:
    template <class data_t>
    using Vars = typename MatterCCZ4RHS<matter_t>::template Vars<data_t>;

    /// derive record name
    static inline const std::string name = "EMTensor";

    /// Variable names
    static inline amrex::Vector<std::string> var_names = {"rho"};

    /// 3 components for the momentum density: Si_11, Si_22, Si_33
    /// 6 components for the stress energy tensor (symmetric): Sij_11, Sij_12
    /// etc.
    static inline const amrex::Vector<std::string> extra_var_names = {
        "Sij_11, Sij_22, Sij_33, Sij_11, Sij_12, Sij_13, Sij_22, Sij_23, "
        "Sij_33"};

    //! Constructor
    EMTensor(const matter_t a_matter, const double dx, const int a_c_rho = -1,
             const Interval a_c_Si  = Interval(),
             const Interval a_c_Sij = Interval());

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute(int i, int j, int k, const amrex::Array4<data_t> &out_mf,
            const amrex::Array4<const data_t> &in_mf) const;

    // Set do_all_components to true to calculate the momentum density
    // and stress energy tensors as well
    AMREX_FORCE_INLINE static void set_up(int a_state_index,
                                          bool do_all_components = false);

    AMREX_FORCE_INLINE static void
    compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
               const amrex::MultiFab &src_mf, const amrex::Geometry &geomdata,
               amrex::Real /*time*/, const int * /*bcrec*/, int /*level*/);

  protected:
    matter_t m_matter;
    FourthOrderDerivatives m_deriv;
    int m_c_rho;      // var enum for the energy density
    Interval m_c_Si;  // Interval of var enums for the momentum density
    Interval m_c_Sij; // Interval of var enums for the spatial
                      // stress-energy density
};

#include "EMTensor.impl.hpp"

#endif /* EMTENSOR_HPP */
