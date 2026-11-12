/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TEUKOLSKYWAVEINITIALDATA_HPP_
#define TEUKOLSKYWAVEINITIALDATA_HPP_

#include "Coordinates.hpp"
#include "EppleyPacket.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

#include <array>
#include <string>

// See the docs for further details on these initial conditions
class TeukolskyWaveInitialData
{
  public:
    //! The parameters of the Teukolsky wave initial data
    struct params_t
    {
        EppleyPacket::params_t packet; //!< parameters of the wave packet
        std::array<amrex::Real, AMREX_SPACEDIM>
            center{}; //!< center of the wave packet

        static void check_params();
        void fill_params();
    };

    //! Reads the parameters, and constructs the wave packet from them
    AMREX_FORCE_INLINE
    explicit TeukolskyWaveInitialData(amrex::Real a_dx);

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const;
    // NOLINTEND(bugprone-easily-swappable-parameters)

  protected:
    amrex::Real m_dx;
    std::array<amrex::Real, AMREX_SPACEDIM> m_center;
    EppleyPacket m_eppley_packet; //!< The Eppley packet object

  private:
    AMREX_FORCE_INLINE
    TeukolskyWaveInitialData(amrex::Real a_dx, const params_t &a_params);

    //! Returns the parameters read from the parameter file
    static params_t read_params();

    //! Converts the parity and azimuthal number m to the type of the packet,
    //! aborting if that combination is not implemented
    static EppleyPacket::Type get_packet_type(const std::string &a_parity,
                                              int a_azimuthal_number);
};

#include "TeukolskyWaveInitialData.impl.hpp"

#endif /* TEUKOLSKYWAVEINITIALDATA_HPP_ */
