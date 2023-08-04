/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DEFAULTLEVELFACTORY_HPP_
#define DEFAULTLEVELFACTORY_HPP_

#include <AMReX_LevelBld.H>

template <class level_t> class DefaultLevelFactory : public amrex::LevelBld
{
  public:
    void variableSetUp() override;
    void variableCleanUp() override;
    amrex::AmrLevel *operator()() override;
    amrex::AmrLevel *
    operator()(amrex::Amr &papa, int lev, const amrex::Geometry &level_geom,
               const amrex::BoxArray &box_array,
               const amrex::DistributionMapping &distribution_mapping,
               amrex::Real time) override;
};

template <class level_t> void DefaultLevelFactory<level_t>::variableSetUp()
{
    level_t::variableSetUp();
}

template <class level_t> void DefaultLevelFactory<level_t>::variableCleanUp()
{
    level_t::variableCleanUp();
}

template <class level_t>
amrex::AmrLevel *DefaultLevelFactory<level_t>::operator()()
{
    // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
    return new level_t;
}

template <class level_t>
amrex::AmrLevel *DefaultLevelFactory<level_t>::operator()(
    amrex::Amr &papa, int lev, const amrex::Geometry &level_geom,
    const amrex::BoxArray &box_array,
    const amrex::DistributionMapping &distribution_mapping, amrex::Real time)
{
    // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
    return new level_t(papa, lev, level_geom, box_array, distribution_mapping,
                       time);
}

#endif /* DEFAULTLEVELFACTORY_HPP_ */
