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

// AMReX Includes
#include <AMReX_MultiFab.H>

//! Calculates the EM tensor and then saves the ones specified in the
//! constructor on the grid

//! Parts of the EMTensor to store
enum class EMTensorOptions
{
    justEnergyDensity,          //! just the energy density
    energyAndMomentumDensities, //! both the energy density and the momentum
                                //! density
    allDensities, //! the energy density, momentum density and the stress tensor
};

//! Calculates the EM tensor and then saves the parts depending on the
//! em_tensor_options_t template argument
template <class matter_t, enum EMTensorOptions em_tensor_options> class EMTensor
{
  public:
    using Vars = typename matter_t::Vars;

    /// derive record name
    static inline const std::string name = "EMTensor";

    /// Variable names
    static amrex::Vector<std::string> var_names();

    //! Constructor
    EMTensor(double dx, int a_dcomp);

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &emtensor_out,
               const amrex::Array4<const amrex::Real> &state) const;
    // NOLINTEND(bugprone-easily-swappable-parameters)

    // Set do_all_components to true to calculate the momentum density
    // and stress energy tensors as well
    AMREX_FORCE_INLINE static void set_up(int a_state_index);

    AMREX_FORCE_INLINE static void
    compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
               const amrex::MultiFab &state_mf, const amrex::Geometry &geomdata,
               amrex::Real /*time*/, const int * /*bcrec*/, int /*level*/);

  protected:
    matter_t m_matter;
    FourthOrderDerivatives m_deriv;
    int m_dcomp; //!< which component in the MultiFab to store the first
                 //!< variable
};

#include "EMTensor.impl.hpp"

#endif /* EMTENSOR_HPP */
