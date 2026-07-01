#ifndef TRAJECTORY_EVALUATOR_HPP_
#define TRAJECTORY_EVALUATOR_HPP_

#include "RLLumpState.hpp"    // RL_MAX_LUMPS
#include "TrajectoryParams.hpp"

#include <algorithm>
#include <cmath>

//! Pure function object: evaluates per-lump independent trajectories at time t.
//!
//! Each lump orbits on its own tilted plane with its own radius, angular
//! velocity and initial phase.  The shared breathing and z-oscillation are
//! superimposed on all lumps.
//!
//! SOLID: Single Responsibility — ONLY trajectory math.
struct TrajectoryEvaluator
{
    //! Evaluate independent trajectories for all active lumps at time ``t``.
    //! Writes Cartesian centres into ``out_x/y/z`` and per-lump pump amplitude
    //! into ``out_amp``.
    //!
    //! Called once per coarse timestep on the CPU.
    static void evaluate(
        double t, const TrajectoryParams &traj,
        double *out_x, double *out_y, double *out_z, double *out_amp)
    {
        const int n = std::min(traj.num_lumps, static_cast<int>(RL_MAX_LUMPS));

        // Shared: radial breathing (added to each lump's R0).
        const double breath_dr =
            traj.A_breath * std::sin(traj.omega_breath * t);

        // Shared: confined axial oscillation.
        const double z_center = traj.z_amp * std::sin(traj.omega_z * t);

        for (int k = 0; k < n; ++k)
        {
            const auto &lk = traj.lumps[k];

            // Per-lump orbital radius with shared breathing and radial drift.
            // Spiral orbits: R0 is the radius at t=0, v_rad is the signed
            // radial drift speed.  Keeps circles as the v_rad=0 special case.
            // Floor at TRAJECTORY_R_MIN so strong inward drift cannot cross r=0.
            const double r_k = std::max(
                TRAJECTORY_R_MIN, lk.R0 + lk.v_rad * t + breath_dr);

            // Per-lump azimuthal angle: phi(t) = phase0 + omega_rot * t
            const double phi = lk.phase0 + lk.omega_rot * t;

            // Position in the lump's own orbital plane (before tilt).
            const double x_orb = r_k * std::cos(phi);
            const double y_orb = r_k * std::sin(phi);
            // z_orb = 0 (orbit lives in a plane).

            // Rotate the orbital plane into lab frame.
            // The orbital-plane normal is defined by (tilt_theta, tilt_phi):
            //   1. Rotate about z by tilt_phi  (picks the tilt direction)
            //   2. Rotate about the new y' by tilt_theta  (tilts the plane)
            //
            // Combined rotation R_z(tilt_phi) * R_y(tilt_theta):
            //   [ cp*ct  -sp   cp*st ] [ x_orb ]
            //   [ sp*ct   cp   sp*st ] [ y_orb ]
            //   [ -st     0    ct    ] [   0   ]
            //
            // where ct=cos(tilt_theta), st=sin(tilt_theta),
            //       cp=cos(tilt_phi),   sp=sin(tilt_phi).
            const double ct = std::cos(lk.tilt_theta);
            const double st = std::sin(lk.tilt_theta);
            const double cp = std::cos(lk.tilt_phi);
            const double sp = std::sin(lk.tilt_phi);

            out_x[k] = cp * ct * x_orb - sp * y_orb;
            out_y[k] = sp * ct * x_orb + cp * y_orb;
            out_z[k] = -st * x_orb + z_center;

            out_amp[k] = std::max(0.0, lk.well_depth);
        }
    }
};

#endif /* TRAJECTORY_EVALUATOR_HPP_ */
