/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef KERRBHINITIALDATA_HPP_
#define KERRBHINITIALDATA_HPP_


#include "CoordinateTransformations.hpp"
#include "Coordinates.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total number of components
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include <AMReX_GpuQualifiers.H>
#include <array>

//! Class which computes the Kerr initial conditions per arXiv 1401.1548
class KerrBHInitialData
{
  public:
    //! Struct for the params of the Kerr BHInitialData
    struct params_t
    {
        amrex::Real mass{}; //!<< The mass of the Kerr BH
        std::array<amrex::Real, AMREX_SPACEDIM> center{}; //!< The center of the Kerr BH
        amrex::Real spin{}; //!< The spin param a = J/M, so 0 <= |a| <= M
        std::array<amrex::Real, AMREX_SPACEDIM> spin_direction = {
            0., 0., 1.}; // default to 'z' axis;
        
        static void check_params();
        inline void fill_params();
    };

  protected:
    amrex::Real m_dx;
    params_t m_params;
    Tensor::Rank2 m_R;

  public:
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_FORCE_INLINE
    KerrBHInitialData(amrex::Real a_dx);
    // NOLINTEND(bugprone-easily-swappable-parameters)

    AMREX_GPU_DEVICE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const;

  protected:
    //! Function which computes the components of the metric in spherical coords
    AMREX_FORCE_INLINE AMREX_GPU_DEVICE void 
    compute_kerr(Tensor::Sym12Rank2 &spherical_g, //!<< The spatial metric in spherical coords
                 Tensor::Sym12Rank2 &spherical_K, //!<< The extrinsic curvature in spherical coords
                 Tensor::Rank1 &spherical_shift, //!<< The spherical components of the shift
                 amrex::Real &kerr_lapse, //!<< The lapse for the kerr solution
                 const Tensor::Rank1 &coords) const;

    //! Helper function to unpack a Sym12Rank2 into a full 3x3 Rank2
    [[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE Tensor::Rank2
    to_rank2(const Tensor::Sym12Rank2& sym) const;
};

#include "KerrBHInitialData.impl.hpp"

#endif /* KERRBHINITIALDATA_HPP_ */
