/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef EMTENSOR_HPP
#define EMTENSOR_HPP

#include "CCZ4Geometry.hpp"
#include "CCZ4RHSWithMatter.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

// AMReX Includes
#include <AMReX_MultiFab.H>

//! Calculates the EM tensor and then saves the ones specified in the
//! constructor on the grid
template <class matter_t> class EMTensor
{
  public:
    template <class data_t>
    using Vars = typename CCZ4RHSWithMatter<matter_t>::template Vars<data_t>;

    /// derive record name
    static inline const std::string name = "EMTensor";

    /// Variable names
    static inline amrex::Vector<std::string> var_names = {"rho"};

    /// 3 components for the momentum density: Si_11, Si_22, Si_33
    /// 6 components for the stress energy tensor (symmetric): Sij_11, Sij_12
    /// etc.
    static inline const amrex::Vector<std::string> extra_var_names = {
        "j_1, j_2, j_3, S_11, S_12, S_13, S_22, S_23, S_33"};

    //! Constructor
    EMTensor(const double dx, const int a_c_rho = -1,
             const Interval a_c_j = Interval(),
             const Interval a_c_S = Interval());

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute(int i, int j, int k, const amrex::Array4<data_t> &out_arrays,
            const amrex::Array4<const data_t> &in_arrays) const;

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
    int m_c_rho;    // var enum for the energy density
    Interval m_c_j; // Interval of var enums for the momentum density
    Interval m_c_S; // Interval of var enums for the spatial
                    // stress-energy density
};

#include "EMTensor.impl.hpp"

#endif /* EMTENSOR_HPP */
