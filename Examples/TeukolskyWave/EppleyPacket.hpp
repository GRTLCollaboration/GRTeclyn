#ifndef EPPLEYPACKET_HPP_
#define EPPLEYPACKET_HPP_

#include "cmath"
#include <AMReX_REAL.H>
#include <vector>

struct EppleyPacketDerivs
{
    amrex::Real F0, F1, F2, F3, F4;
};

struct EppleyPacketMetricComponents
{
    amrex::Real gxx, gxy, gxz, gyy, gyz, gzz;
};

//! Base EppleyPacket class
class EppleyPacket
{
  public:
    struct params_t
    {
        amrex::Real sigma{};         //!< width of the packet
        amrex::Real amplitude{};     //!< amplitude of the packet
        amrex::Real radial_offset{}; // offset for radial coordinate to not
                                     // center the wave on the center
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
    }

    params_t m_params;

    //! F function and its first four derivatives, where x = r \pm t
    EppleyPacketDerivs get_F_derivs(amrex::Real x) const;

    EppleyPacket() { m_params.fill_params(); };

    EppleyPacketMetricComponents
    get_metric_components(amrex::Real x, amrex::Real y, amrex::Real z) const;
};

struct EvenEppleyPacketCoefficients
{
    amrex::Real A, B, C;
};
//! Base class for the even parity Eppley packets
class EvenEppleyPacket : EppleyPacket
{
  public:
    EvenEppleyPacketCoefficients get_ABC(amrex::Real r) const;

    EvenEppleyPacket() : EppleyPacket() {};
};

struct OddEppleyPacketCoefficients
{
    amrex::Real K, L;
};
//! Base class for the odd parity Eppley packets
class OddEppleyPacket : EppleyPacket
{
  public:
    OddEppleyPacketCoefficients get_KL(amrex::Real r) const;

    OddEppleyPacket() : EppleyPacket() {};
};

//! Specific Eppley Packet classes
//! m = 0 and m = 2 are the only ones implemented so far for the even parity
//! case
class EvenEppleyPacketM0 : public EvenEppleyPacket
{
  public:
    EvenEppleyPacketM0() : EvenEppleyPacket() {}

    EppleyPacketMetricComponents
    get_metric_components(amrex::Real x, amrex::Real y,
                          amrex::Real z) const override;
};

class EvenEppleyPacketM2 : public EvenEppleyPacket
{
  public:
    EvenEppleyPacketM2() : EvenEppleyPacket() {}

    EppleyPacketMetricComponents
    get_metric_components(amrex::Real x, amrex::Real y,
                          amrex::Real z) const override;
};

//! m = 2 class for the odd parity case
class OddEppleyPacketM2 : public OddEppleyPacket
{
  public:
    OddEppleyPacketM2() : OddEppleyPacket() {}

    EppleyPacketMetricComponents
    get_metric_components(amrex::Real x, amrex::Real y,
                          amrex::Real z) const override;
};

#include "EppleyPacket.impl.hpp"

#endif /* EPPLEYPACKET_HPP_ */
