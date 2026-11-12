/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(TEUKOLSKYWAVEINITIALDATA_HPP_)
#error "This file should only be included through TeukolskyWaveInitialData.hpp"
#endif

#ifndef TEUKOLSKYWAVEINITIALDATA_IMPL_HPP_
#define TEUKOLSKYWAVEINITIALDATA_IMPL_HPP_

#include "GRParmParse.hpp"
#include "TensorAlgebra.hpp"
#include "TeukolskyWaveInitialData.hpp"

#include <cmath>

inline EppleyPacket::Type
TeukolskyWaveInitialData::get_packet_type(const std::string &a_parity,
                                          int a_azimuthal_number)
{
    if (a_parity == "even" && a_azimuthal_number == 0)
    {
        return EppleyPacket::Type::even_m0;
    }
    if (a_parity == "even" && a_azimuthal_number == 2)
    {
        return EppleyPacket::Type::even_m2;
    }
    if (a_parity == "odd" && a_azimuthal_number == 2)
    {
        return EppleyPacket::Type::odd_m2;
    }
    amrex::Abort("Combination of parity and azimuthal number not implemented."
                 " Must be (even, 0), (even, 2) or (odd, 2).");
    return EppleyPacket::Type::even_m0; // never reached
}

inline void TeukolskyWaveInitialData::params_t::check_params()
{
    GRParmParse tw_pp("teukolsky_wave");

    // Check that the particular choice of parity and m is implemented
    int azimuthal_number{0};
    std::string parity{"even"};
    tw_pp.queryAdd("azimuthal_number", azimuthal_number);
    tw_pp.queryAdd("parity", parity);
    get_packet_type(parity, azimuthal_number);

    // These are parameters for the seed function
    amrex::Real sigma{};
    tw_pp.queryAdd("sigma", sigma);
    if (sigma <= 0.0)
    {
        tw_pp.error("sigma", "must be > 0");
    }
    amrex::Real amplitude{};
    tw_pp.queryAdd("amplitude", amplitude);
    if (amplitude <= 0.0)
    {
        tw_pp.error("amplitude", "must be > 0");
    }
    amrex::Real radial_offset{};
    tw_pp.queryAdd("radial_offset", radial_offset);
    if (radial_offset <= 0.0)
    {
        tw_pp.error("radial_offset", "must be > 0");
    }
    amrex::Real regularize_r{};
    tw_pp.queryAdd("regularize_r", regularize_r);
    if (regularize_r <= 0.0)
    {
        tw_pp.error("regularize_r", "must be > 0");
    }
}

inline void TeukolskyWaveInitialData::params_t::fill_params()
{
    GRParmParse tw_pp("teukolsky_wave");
    int azimuthal_number{};
    std::string parity{};
    tw_pp.get("azimuthal_number", azimuthal_number);
    tw_pp.get("parity", parity);
    packet.type = get_packet_type(parity, azimuthal_number);
    tw_pp.get("sigma", packet.sigma);
    tw_pp.get("amplitude", packet.amplitude);
    tw_pp.get("radial_offset", packet.radial_offset);
    tw_pp.get("regularize_r", packet.regularize_r);

    GRParmParse geometry_pp("geometry");
    geometry_pp.get("center", center);
}

inline TeukolskyWaveInitialData::params_t
TeukolskyWaveInitialData::read_params()
{
    params_t params;
    params.fill_params();
    return params;
}

AMREX_FORCE_INLINE
TeukolskyWaveInitialData::TeukolskyWaveInitialData(amrex::Real a_dx)
    : TeukolskyWaveInitialData(a_dx, read_params())
{
}

AMREX_FORCE_INLINE
TeukolskyWaveInitialData::TeukolskyWaveInitialData(amrex::Real a_dx,
                                                   const params_t &a_params)
    : m_dx(a_dx), m_center(a_params.center), m_eppley_packet(a_params.packet)
{
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void TeukolskyWaveInitialData::operator()(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &state) const
{
    const amrex::CellData<amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    const Coordinates coords(amrex::IntVect(ix, iy, iz), m_dx, m_center);

    const Tensor::Sym12Rank2 g =
        m_eppley_packet.get_metric(coords.x, coords.y, coords.z);

    // Define the conformal factor
    const amrex::Real chi =
        std::pow(TensorAlgebra::compute_determinant_sym(g), -1. / 3.);

    FOR2_SYM(i, j)
    {
        state_cell_data[sym_var_idx(c_h11, i, j)] = chi * g(i, j);
    }

    state_cell_data[c_chi]   = chi;
    state_cell_data[c_lapse] = 1.;

    // This will need to be set by the GammaCalculator later, since the metric
    // is not conformally flat, but we will zero it here for now.
    FOR (i)
    {
        state_cell_data[c_Gamma1 + i] = 0.0_rt;
    }

    // TeukolskyWaveLevel zero-initializes every state component before calling
    // this operator, so K, A_ij, Theta, shift and B remain zero here. This
    // is consistent with the time-symmetric nature of the data.
}

#endif /* TEUKOLSKYWAVEINITIALDATA_IMPL_HPP_ */
