/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TEUKOLSKYWAVEINITIALDATA_HPP_
#define TEUKOLSKYWAVEINITIALDATA_HPP_

#include "Coordinates.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"
#include "EppleyPacket.hpp"

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

#include <array>
#include <cstddef>

// See the docs for further details on these initial conditions
class TeukolskyWaveInitialData
{
  public:
    struct params_t
    {

        static int magnetic{};       //!< magnetic quantum number
        static std::string parity{}; //!< parity of the packet. 0 for even, 1 for odd
	static amrex::Real sigma{};  //!< width of the packet
        static amrex::Real amplitude{};     //!< amplitude of the packet
        static amrex::Real radial_offset{}; // offset for radial coordinate to not
                                     // center the wave on the center
        static amrex::Real regularize_r{};  // small regularization parameter for the
                                     // radial coordinate
        static void check_params()
        {
            // These are parameters specfic to the Teukolsky wave example
            GRParmParse tw_pp("teukolsky_wave");
            tw_pp.queryAdd("sigma", sigma);
            if (sigma <= 0.0)
            {
                tw_pp.error("sigma", "must be > 0");
            }
            tw_pp.queryAdd("amplitude", amplitude);
            if (amplitude <= 0.0)
            {
                tw_pp.error("amplitude", "must be > 0");
            }
            tw_pp.queryAdd("radial_offset", radial_offset);
            if (radial_offset < 0.0)
            {
                tw_pp.error("radial_offset", "must be >= 0");
            }
            tw_pp.queryAdd("regularize_r", regularize_r);
            if (regularize_r <= 0.0)
            {
                tw_pp.error("regularize_r", "must be > 0");
            }
        }

        void fill_params()
        {
            GRParmParse tw_pp("teukolsky_wave");
            tw_pp.get("sigma", sigma);
            tw_pp.get("amplitude", amplitude);
            tw_pp.get("radial_offset", radial_offset);
            tw_pp.get("regularize_r", regularize_r);
        }
    };

    AMREX_FORCE_INLINE
    explicit TeukolskyWaveInitialData(amrex::Real a_dx);

    AMREX_FORCE_INLINE
    void initialize_eppley_packet(int magnetic, std::string parity);

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const;
    // NOLINTEND(bugprone-easily-swappable-parameters)

  protected:
    amrex::Real m_dx;
    std::array<amrex::Real, AMREX_SPACEDIM> m_center{};
    params_t m_params;
    EppleyPacket m_eppley_packet; //!< The Eppley packet object
};

#include "TeukolskyWaveInitialData.impl.hpp"

#endif /* TEUKOLSKYWAVEINITIALDATA_HPP_ */
