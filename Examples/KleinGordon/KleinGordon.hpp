#ifndef KLEINGORDON_HPP_
#define KLEINGORDON_HPP_

#include "GRAMR.hpp"

class KleinGordon : public GRAMR
{
  public:

    KleinGordon(amrex::LevelBld *a_levelbld) : GRAMR(a_levelbld) {}

    // If we need to override any virtual functions in amrex::Amr, we can do
    // it here.
    virtual ~KleinGordon()
    {
        amrex::MultiFab const &state_new =
            this->getLevel(0).get_new_data(State_Type);
        amrex::Print()
            << "At the end of simulation, the min and max of the wave are "
            << state_new.min(0) << " and " << state_new.max(0) << "\n\n";
    };
};

#endif /* KLEINGORDON_HPP */
