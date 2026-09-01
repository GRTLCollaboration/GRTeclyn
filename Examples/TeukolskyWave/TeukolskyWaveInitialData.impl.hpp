/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(TEUKOLSKYWAVEINITIALDATA_HPP_)
#error "This file should only be included through TeukolskyWaveInitialData.hpp"
#endif

#ifndef TEUKOLSKYWAVEINITIALDATA_IMPL_HPP_
#define TEUKOLSKYWAVEINITIALDATA_IMPL_HPP_

#include "EppleyPacket.hpp"
#include "GRParmParse.hpp"
#include "TensorAlgebra.hpp"
#include "TeukolskyWaveInitialData.hpp"

#include <cmath>

AMREX_FORCE_INLINE
TeukolskyWaveInitialData::TeukolskyWaveInitialData(amrex::Real a_dx)
    : m_dx(a_dx)
{
    m_params.fill_params();
    // Initialize the center of the geometry for the wave
    GRParmParse geometry_pp("geometry");
    geometry_pp.get("center", m_center);
}

void TeukolskyWaveInitialData::initialize_eppley_packet(int magnetic,
                                                        std::string parity)
{
    m_params.magnetic = magnetic;
    m_params.parity   = parity;
    if (parity == "even" && magnetic == 0)
    {
        m_eppley_packet = EvenEppleyPacketM0();
    }
    else if (parity == "even" && magnetic == 2)
    {
        m_eppley_packet = EvenEppleyPacketM2();
    }
    else if (parity == "odd" && magnetic == 2)
    {
        m_eppley_packet = OddEppleyPacketM2();
    }
    else
    {
        amrex::Abort("Invalid combination of parity and magnetic quantum "
                     "number. Must be (even, 0), (even, 2) or (odd, 2).");
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void TeukolskyWaveInitialData::operator()(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &state) const
{
    const amrex::CellData<amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    const Coordinates coords(amrex::IntVect(ix, iy, iz), m_dx, m_center);

    const EppleyPacketMetricComponents metric =
        m_eppley_packet.get_metric_components(coords.x, coords.y, coords.z);

    const amrex::Real g[3][3] = {
        {metric.gxx, metric.gxy, metric.gxz},
        {metric.gxy, metric.gyy, metric.gyz},
        {metric.gxz, metric.gyz, metric.gzz}
    };

    // Define the conformal factor
    amrex::Real det_g = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1]) -
                        g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0]) +
                        g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
    amrex::Real chi = pow(det_g, -1. / 3.);

    FOR2_SYM(i, j)
    {
        state_cell_data[sym_var_idx(c_h11, i, j)] = chi * g[i][j];
    }

    state_cell_data[c_chi]   = chi;
    state_cell_data[c_lapse] = 1.;

    // BELOW ????
    // // This will need to be set by the GammaCalculator, since the metric
    // // is not conformally flat, but we will zero it here for now.
    // FOR (i)
    // {
    //     state_cell_data[c_Gamma1 + i] = 0.0_rt;
    // }

    // TeukolskyWaveLevel zero-initializes every state component before calling
    // this operator, so K, A_ij, Theta, shift and B remain zero here.
}

#endif /* TeukolskyWaveInitialData_IMPL_HPP_ */
