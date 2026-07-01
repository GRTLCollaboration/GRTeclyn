#ifndef TRAJECTORY_PARAMS_HPP_
#define TRAJECTORY_PARAMS_HPP_

#include "RLLumpState.hpp" // RL_MAX_LUMPS

//! Minimum orbital radius used by TrajectoryEvaluator (geometric units).
//! Prevents inward spirals from crossing the origin and inverting the orbit.
inline constexpr double TRAJECTORY_R_MIN = 0.1;

//! Per-lump trajectory: each lump orbits independently on its own tilted
//! plane with its own radius, angular velocity, initial phase, pump
//! amplitude, and optional radial drift.  Counter-rotating lumps (opposite
//! omega_rot signs) create shear and frame dragging; radial drift creates
//! spiral orbits and transient close approaches — the primary FTL mechanisms.
//!
//! SOLID: Single Responsibility — this struct is ONLY the per-lump data.
struct PerLumpTrajectory
{
    double R0{5.0};           //!< orbital radius at t=0
    double omega_rot{0.0};    //!< angular velocity  (sign = direction)
    double phase0{0.0};       //!< initial azimuthal angle [rad]
    double tilt_theta{0.0};   //!< orbital-plane tilt from z-axis [rad]
    double tilt_phi{0.0};     //!< orbital-plane azimuth [rad]
    double well_depth{0.05};  //!< pump amplitude for this lump
    double v_rad{0.0};        //!< radial drift speed (signed: inward < 0, outward > 0)
};

//! Shared trajectory parameters + array of per-lump orbits.
//!
//! The separation is physically motivated: per-lump parameters define
//! independent orbit geometry (where each matter lump circles), while
//! shared parameters control collective phenomena (breathing, z-oscillation,
//! pump width) that would be redundant to duplicate per lump.
struct TrajectoryParams
{
    // --- Per-lump orbits ---
    int num_lumps{0};
    PerLumpTrajectory lumps[RL_MAX_LUMPS]{};

    // --- Shared: radial breathing (applied to ALL lump radii) ---
    double A_breath{0.0};      //!< pulsation amplitude (added to each R0_k)
    double omega_breath{0.0};  //!< pulsation angular frequency

    // --- Shared: confined axial oscillation ---
    //! Replaces the old z_drift*t (which diverged).
    //! z_center(t) = z_amp * sin(omega_z * t)
    //! Always bounded: |z_center| <= z_amp.
    double z_amp{0.0};         //!< axial oscillation amplitude
    double omega_z{0.0};       //!< axial oscillation frequency

    // --- Shared: pump width ---
    double well_width{1.5};    //!< Gaussian spotlight sigma (all lumps)
};

#endif /* TRAJECTORY_PARAMS_HPP_ */
