#include "KleinGordon.hpp"
#include "KleinGordonLevel.hpp"

using namespace amrex;

KleinGordon::~KleinGordon ()
{
    MultiFab const& S = this->getLevel(0).get_new_data(State_Type);
    amrex::Print() << "At the end of simulation, the min and max of the wave are "
                   << S.min(0) << " and " << S.max(0) << "\n\n";
}
