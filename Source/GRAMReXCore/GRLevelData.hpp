/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRLEVELDATA_HPP_
#define GRLEVELDATA_HPP_

#include <AMReX_MultiFab.H>

class GRLevelData : public amrex::MultiFab
{
  public:
    GRLevelData();

    void setVal(const double a_val);

    void setVal(const double a_val, const int a_comp);

//xxxxx    void setVal(const double a_val, const Interval a_comps);

    // loop only goes over a_disjoint_box_layout
//xxxxx    void plus(const GRLevelData &a_src, const double a_scale,
//              const DisjointBoxLayout &a_disjoint_box_layout);
};

#endif /* GRLEVELDATA_HPP_ */
