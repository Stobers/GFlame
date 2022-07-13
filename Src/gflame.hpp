#ifndef _gflame_hpp_
#define _gflame_hpp_

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

#define Pi 3.1415926535897932384626433

namespace Gflame{
    //
    // variables 
    //
    static int max_grid_size = 16;
    static int nsteps = -1;
    static amrex::Real end_time = -1;
    static amrex::Vector<int> ncells, is_periodic;
    static amrex::Vector<amrex::Real> prob_lo, prob_hi;
    static amrex::DistributionMapping dm;
    static amrex::BoxArray ba;
    static amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dx;
    static amrex::MultiFab feild_new, feild;
    static amrex::Geometry geom;
    static const unsigned int feild_ncomp = 2;
    static const unsigned int feild_ngrow = 0;
    static amrex::Real flameloc = 0.5;
    static amrex::Real flamepert = 0;
   
    //
    // functions
    //
    void init();
    void set_bc();
    void init_gfeild();
}
#endif
