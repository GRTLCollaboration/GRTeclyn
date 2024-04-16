#ifndef KLEINGORDON_HPP_
#define KLEINGORDON_HPP_

#include "GRAMR.hpp"

class KleinGordon
    : public GRAMR
{
public:

  KleinGordon(amrex::LevelBld *a_levelbld): GRAMR(a_levelbld){}

  //    virtual ~KleinGordon ();

    // If we need to override any virtual functions in amrex::Amr, we can do
    // it here.


};

#endif   /* KLEINGORDON_HPP */
