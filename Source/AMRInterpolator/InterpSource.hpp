/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INTERPSOURCE_H_
#define INTERPSOURCE_H_

// Other inclues
#include "VariableType.hpp"

#include <AMReX_MultiFab.H>

#include <array>

// Abstract base class to get the FABs out of an AMRLevel
class InterpSource
{
  public:
    [[nodiscard]] virtual const amrex::MultiFab &getLevelData(
        const VariableType var_type = VariableType::evolution) const = 0;
    [[nodiscard]] virtual bool
    contains(const std::array<double, AMREX_SPACEDIM> &point) const = 0;

    virtual ~InterpSource() = default;
};

#endif /* INTERPSOURCE_H_ */
