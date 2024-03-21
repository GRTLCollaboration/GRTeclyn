/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef DEBUGGINGTOOLS_HPP_
#define DEBUGGINGTOOLS_HPP_

#ifdef EQUATION_DEBUG_MODE
#include <AMReX_IntVect.H>
#endif

// Other includes
#include <string.h>

/// This file contains a collection of helpful #defines and other definitions
/// that are hepful for debugging.

// Unfortunately, most of the functionality can only be achieved with macros
// (e.g. including the variable and filename).

#define __FILENAME__                                                           \
    (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#define DEBUG_SHOW(VAR) amrex::Print() << #VAR << ": " << VAR << " "
#define DEBUG_FILE amrex::Print() << __FILENAME__ << ": "
#define DEBUG_END amrex::Print() << std::endl
#define DEBUG_DOUBLE_PRECISION amrex::Print() << std::setprecision(16)

/// The macros DEBUG_OUT make debugging quicker and allow easy printing of a
/// variable.
#define DEBUG_OUT(VAR)                                                         \
    DEBUG_FILE;                                                                \
    DEBUG_SHOW(VAR);                                                           \
    DEBUG_END
#define DEBUG_OUT2(VAR1, VAR2)                                                 \
    DEBUG_FILE;                                                                \
    DEBUG_SHOW(VAR1);                                                          \
    DEBUG_SHOW(VAR2);                                                          \
    DEBUG_END
#define DEBUG_OUT3(VAR1, VAR2, VAR3)                                           \
    DEBUG_FILE;                                                                \
    DEBUG_SHOW(VAR1);                                                          \
    DEBUG_SHOW(VAR2);                                                          \
    DEBUG_SHOW(VAR3);                                                          \
    DEBUG_END
#define DEBUG_OUT4(VAR1, VAR2, VAR3, VAR4)                                     \
    DEBUG_FILE;                                                                \
    DEBUG_SHOW(VAR1);                                                          \
    DEBUG_SHOW(VAR2);                                                          \
    DEBUG_SHOW(VAR3);                                                          \
    DEBUG_SHOW(VAR4);                                                          \
    DEBUG_END

//
#ifdef EQUATION_DEBUG_MODE
#define DEBUG_HEADER                                                           \
    amrex::Print() << "Debug output in " << __FILENAME__                       \
                   << " at: " << s_current_integer_coords << "." << std::endl
static amrex::IntVect s_current_integer_coords;
namespace EquationDebugging
{
inline void check_no_omp()
{
#ifdef _OPENMP
    if (omp_get_max_threads() > 1)
        amrex::Abort("Equation debug mode can  only be used with one thread.");
#endif
}

inline void
set_global_cell_coordinates(const amrex::IntVect current_integer_coords)
{
    check_no_omp();
    s_current_integer_coords = current_integer_coords;
}
} // namespace EquationDebugging
#endif

#endif /* DEBUGGINGTOOLS_HPP_ */
