#ifndef RL_LUMP_STATE_HPP_
#define RL_LUMP_STATE_HPP_

#include "GRTresnaScalarLayout.hpp" // GRTRESNA_MAX_INDEPENDENT_SCALARS

//! Live, per-lump 3-D state measured by the RL lump tracker each action
//! interval.  One of these per lump is exposed to the agent (observation) and
//! its centre drives the matching pump "spotlight" (the spotlight follows the
//! lump).  All quantities are physical (code units), measured on the finest
//! level over a search ball around the lump's previous tracked centre.
struct RLLumpState
{
    double x{0.0};         //!< mass-weighted centroid x (physical coords)
    double y{0.0};         //!< mass-weighted centroid y
    double z{0.0};         //!< mass-weighted centroid z
    double size{0.0};      //!< RMS radius about the centroid (spread / dispersal)
    double mass{0.0};      //!< integral of the weight field (|Phi|^2) in the ball
    double peak{0.0};      //!< max weight in the ball (compactness proxy)
    double min_lapse{1.0}; //!< min lapse in the ball (per-lump collapse warning)
    double min_chi{1.0};   //!< min chi in the ball (per-lump horizon warning)
};

//! Number of doubles each lump contributes to the observation vector.
inline constexpr int RL_LUMP_OBS_STRIDE = 8;

//! Max lumps the tracker / pump can address (shares the matter slot budget).
inline constexpr int RL_MAX_LUMPS = GRTRESNA_MAX_INDEPENDENT_SCALARS;

#endif /* RL_LUMP_STATE_HPP_ */
