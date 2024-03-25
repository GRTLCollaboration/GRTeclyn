#ifndef KLEINGORDON_HPP_
#define KLEINGORDON_HPP_

#include <AMReX_Amr.H>

class KleinGordon
    : public amrex::Amr
{
public:

    using amrex::Amr::Amr;

    virtual ~KleinGordon ();

    // If we need to override any virtual functions in amrex::Amr, we can do
    // it here.


};

#endif
