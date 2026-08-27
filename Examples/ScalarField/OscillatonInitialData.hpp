/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef OSCILLATONINITIALDATA_HPP_
#define OSCILLATONINITIALDATA_HPP_

#include "Coordinates.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

#include <array>
#include <cstddef>

using namespace amrex::literals;

// See the docs for further details on these initial conditions
class OscillatonInitialData
{
  public:
    AMREX_FORCE_INLINE
    explicit OscillatonInitialData(amrex::Real a_dx);

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const;
    // NOLINTEND(bugprone-easily-swappable-parameters)

  protected:
    amrex::Real m_dx;
    std::array<amrex::Real, AMREX_SPACEDIM> m_center{};

    // These parameters come from fitting a known oscillaton solution
    static constexpr amrex::Real geometry_scale = 7.0_rt;
    static constexpr amrex::Real scalar_scale   = 10.0_rt;
    static constexpr amrex::Real scalar_central_value =
        -0.044159713033393284_rt;

    static constexpr int num_compactness_coefficients      = 13;
    static constexpr int num_lapse_correction_coefficients = 11;
    static constexpr int num_scalar_exponent_coefficients  = 13;

    std::array<amrex::Real, num_compactness_coefficients>
        m_compactness_coefficients{
            2.8958846304099733e-01_rt,  -1.3137514414560530e-01_rt,
            -1.6701297363810957e-02_rt, 6.0707390156588582e-03_rt,
            2.3289534559576405e-03_rt,  2.8137187117568101e-04_rt,
            -1.2773355073195884e-04_rt, -1.0740930032253914e-04_rt,
            -5.9245056266891001e-05_rt, -1.1536669472411421e-05_rt,
            -5.4379421632360992e-06_rt, 5.0832422596546345e-06_rt,
            -1.4123119616337690e-06_rt};

    std::array<amrex::Real, num_lapse_correction_coefficients>
        m_lapse_correction_coefficients{
            -1.1599758613808658e-01_rt, 1.4759576543983041e-01_rt,
            -1.5983935412197312e-02_rt, -1.5671775506787750e-02_rt,
            -3.7118770954652227e-04_rt, 1.9275458497198111e-03_rt,
            9.4289642098836898e-04_rt,  2.2891955428106675e-04_rt,
            -3.5029665945347518e-05_rt, -7.1034651695464193e-05_rt,
            -3.3347291268122290e-05_rt};

    std::array<amrex::Real, num_scalar_exponent_coefficients>
        m_scalar_exponent_coefficients{
            3.0843643066252771e00_rt,  2.8024295695301471e-01_rt,
            1.7478099724814108e-01_rt, 1.4662523126960619e-01_rt,
            9.7285887654544018e-02_rt, 6.8961110591586425e-02_rt,
            4.8276770778464026e-02_rt, 3.0232533883913594e-02_rt,
            1.8394483397716436e-02_rt, 9.3598232915395006e-03_rt,
            4.5638233137690540e-03_rt, 1.5222932582614564e-03_rt,
            4.9664880641856293e-04_rt};

    template <std::size_t N>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    evaluate_chebyshev(amrex::Real x,
                       const std::array<amrex::Real, N> &coefficients) const;
};

#include "OscillatonInitialData.impl.hpp"

#endif /* OSCILLATONINITIALDATA_HPP_ */
