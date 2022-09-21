/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ALWAYS_INLINE_HPP_
#define ALWAYS_INLINE_HPP_

#include <AMReX_Extension.H>

#define ALWAYS_INLINE AMREX_FORCE_INLINE

#if 0
#if defined(__GNUC__)
#define ALWAYS_INLINE __attribute__((always_inline)) __inline__
#else
#define ALWAYS_INLINE inline
#endif
#endif

#endif /* ALWAYS_INLINE_HPP_ */
