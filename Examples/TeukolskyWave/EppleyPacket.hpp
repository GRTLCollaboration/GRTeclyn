/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef EPPLEYPACKET_HPP_
#define EPPLEYPACKET_HPP_

#include "Tensor.hpp"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <cmath>

//! Superposition of an ingoing and outgoing Teukolsky wave.
class EppleyPacket
{
  public:
    //! Which (parity, azimuthal number m) combination a packet represents
    enum class Type
    {
        even_m0,
        even_m2,
        odd_m2,
    };

    //! The parameters of the packet (read by TeukolskyWaveInitialData)
    struct params_t
    {
        Type type{Type::even_m0};
        amrex::Real sigma{};         //!< width of the packet
        amrex::Real amplitude{};     //!< amplitude of the packet
        amrex::Real radial_offset{}; //!< offset for the radial coordinate, so
                                     //!< that the wave is not centred at the
                                     //!< origin
        amrex::Real regularize_r{};  //!< small regularization parameter for the
                                     //!< radial coordinate
    };

    explicit EppleyPacket(const params_t &a_params) : m_params(a_params) {}

    //! Get the spatial metric of the Eppley packet at t = 0
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
    get_metric(amrex::Real x, amrex::Real y, amrex::Real z) const;

  private:
    //! Metric coefficients used in an even parity Teukolsky wave
    struct EvenCoefficients
    {
        amrex::Real A, B, C;
    };

    //! Metric coefficients used in an odd parity Teukolsky wave
    struct OddCoefficients
    {
        amrex::Real K, L;
    };

    //! The seed function F and its first four derivatives
    struct SeedDerivs
    {
        amrex::Real F0, F1, F2, F3, F4;
    };

    params_t m_params;

    //! m = 0 even parity case
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
    get_metric_even_m0(amrex::Real x, amrex::Real y, amrex::Real z) const;

    //! m = 2 even parity case
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
    get_metric_even_m2(amrex::Real x, amrex::Real y, amrex::Real z) const;

    //! m = 2 odd parity case
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
    get_metric_odd_m2(amrex::Real x, amrex::Real y, amrex::Real z) const;

    //! Get F and its first four derivatives, where x = r \pm t
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE SeedDerivs
    get_F_derivs(amrex::Real x) const;

    //! Get the metric coefficients A, B, and C for the even parity case
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE EvenCoefficients
    get_ABC(amrex::Real r) const;

    //! Get the metric coefficients K and L for the odd parity case
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE OddCoefficients
    get_KL(amrex::Real r) const;

    //! Numerical regularization of the radial coordinate, to avoid dividing
    //! by zero at the origin (see the docs). It only changes r for r <~
    //! regularize_r.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    regularize_r(amrex::Real r) const
    {
        return r + m_params.regularize_r * std::exp(-r * r /
                                                    (m_params.regularize_r *
                                                     m_params.regularize_r));
    }
};

#include "EppleyPacket.impl.hpp"

#endif /* EPPLEYPACKET_HPP_ */
