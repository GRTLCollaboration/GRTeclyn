#ifndef EPPLEYPACKET_HPP_
#define EPPLEYPACKET_HPP_

#include "GRParmParse.hpp"

#include <AMReX_REAL.H>
#include <cmath>

//! Metric coefficients used in an even parity Teukolsky wave
struct EvenEppleyPacketCoefficients
{
    amrex::Real A, B, C;
};

//! Metric coefficients used in an odd parity Teukolsky wave
struct OddEppleyPacketCoefficients
{
    amrex::Real K, L;
};

//! Derivatives of the seed function F that are used in the metric coefficients
//! of a Teukolsky wave
struct EppleyPacketDerivs
{
    amrex::Real F0, F1, F2, F3, F4;
};

//! Struct to wrap the metric components of the Eppley packet
struct EppleyPacketMetricComponents
{
    amrex::Real gxx, gxy, gxz, gyy, gyz, gzz;
};

//! Which (parity, magnetic quantum number) combination a packet represents
enum class EppleyPacketType
{
    even_m0,
    even_m2,
    odd_m2
};

//! Superposition of an ingoing and outgoing Teukolsky wave.
//! This is captured by value into GPU device lambdas (via
//! TeukolskyWaveInitialData), so it must stay trivially copyable: the choice
//! of (parity, magnetic number) is a plain enum tag.
class EppleyPacket
{
  public:
    struct params_t
    {
        amrex::Real sigma{};         //!< width of the packet
        amrex::Real amplitude{};     //!< amplitude of the packet
        amrex::Real radial_offset{}; // offset for radial coordinate to not
                                     // center the wave at the origin
        amrex::Real regularize_r{};  // small regularization parameter for the
                                     // radial coordinate
        void fill_params()
        {
            GRParmParse tw_pp("teukolsky_wave");
            tw_pp.get("sigma", sigma);
            tw_pp.get("amplitude", amplitude);
            tw_pp.get("radial_offset", radial_offset);
            tw_pp.get("regularize_r", regularize_r);
        }
    };

    params_t m_params;
    EppleyPacketType m_type{EppleyPacketType::even_m0};

    EppleyPacket() { m_params.fill_params(); }

    explicit EppleyPacket(EppleyPacketType a_type) : m_type(a_type)
    {
        m_params.fill_params();
    }

    //! F function and its first four derivatives, where x = r \pm t
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE EppleyPacketDerivs
    get_F_derivs(amrex::Real x) const;

    [[nodiscard]] AMREX_GPU_DEVICE
        AMREX_FORCE_INLINE EvenEppleyPacketCoefficients
        get_ABC(amrex::Real r) const;

    [[nodiscard]] AMREX_GPU_DEVICE
        AMREX_FORCE_INLINE OddEppleyPacketCoefficients
        get_KL(amrex::Real r) const;

    [[nodiscard]] AMREX_GPU_DEVICE
        AMREX_FORCE_INLINE EppleyPacketMetricComponents
        get_metric_components(amrex::Real x, amrex::Real y,
                              amrex::Real z) const;

  private:
    //! m = 0 even parity case
    [[nodiscard]] AMREX_GPU_DEVICE
        AMREX_FORCE_INLINE EppleyPacketMetricComponents
        get_metric_components_even_m0(amrex::Real x, amrex::Real y,
                                      amrex::Real z) const;

    //! m = 2 even parity case
    [[nodiscard]] AMREX_GPU_DEVICE
        AMREX_FORCE_INLINE EppleyPacketMetricComponents
        get_metric_components_even_m2(amrex::Real x, amrex::Real y,
                                      amrex::Real z) const;

    //! m = 2 odd parity case
    [[nodiscard]] AMREX_GPU_DEVICE
        AMREX_FORCE_INLINE EppleyPacketMetricComponents
        get_metric_components_odd_m2(amrex::Real x, amrex::Real y,
                                     amrex::Real z) const;
};

#include "EppleyPacket.impl.hpp"

#endif /* EPPLEYPACKET_HPP_ */
