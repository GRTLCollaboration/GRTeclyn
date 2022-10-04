/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTERPSOURCE_H_
#define INTERPSOURCE_H_

// Other inclues
#include "VariableType.hpp"
#include <array>

#include <AMReX_MultiFab.H>

// Abstrace base class to get the FABs out of an AMRLevel
class InterpSource
{
  public:
    virtual const MultiFab &getLevelData(
        const VariableType var_type = VariableType::evolution) const = 0;
    virtual bool
    contains(const std::array<double, AMREX_SPACEDIM> &point) const = 0;
};

#endif /* INTERPSOURCE_H_ */
