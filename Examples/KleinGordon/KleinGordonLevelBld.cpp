#include "KleinGordonLevel.hpp"
#include <AMReX_LevelBld.H>

using namespace amrex;

class KleinGordonLevelBld
    : public LevelBld
{
    virtual void variableSetUp () override;
    virtual void variableCleanUp () override;
    virtual AmrLevel *operator() () override;
    virtual AmrLevel *operator() (Amr& amr, int lev, const Geometry& gm,
                                  const BoxArray& ba, const DistributionMapping& dm,
                                  Real time) override;
};

KleinGordonLevelBld KleinGordon_bld;

LevelBld*
getLevelBld ()
{
    return &KleinGordon_bld;
}

void
KleinGordonLevelBld::variableSetUp ()
{
    KleinGordonLevel::variableSetUp();
}

void
KleinGordonLevelBld::variableCleanUp ()
{
    KleinGordonLevel::variableCleanUp();
}

AmrLevel*
KleinGordonLevelBld::operator() ()
{
    return new KleinGordonLevel;
}

AmrLevel*
KleinGordonLevelBld::operator() (Amr& amr, int lev, const Geometry& gm,
                          const BoxArray& ba, const DistributionMapping& dm,
                          Real time)
{
    return new KleinGordonLevel(amr, lev, gm, ba, dm, time);
}
