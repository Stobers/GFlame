#include <gflame.hpp>

#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>


namespace Gflame{
    void init()
    {
	amrex::Print() << "... initialising domain \n";
	
	amrex::ParmParse pp_dom("dom");
	pp_dom.getarr("ncells",ncells);
	pp_dom.getarr("prob_lo",prob_lo);
	pp_dom.getarr("prob_hi",prob_hi);
	pp_dom.getarr("is_periodic",is_periodic);
	pp_dom.query("max_grid_size",max_grid_size);

	amrex::ParmParse pp;
	pp.query("nsteps",nsteps);
	pp.query("end_time",end_time);
	if (nsteps < 0 && end_time < 0)
	    amrex::Abort("both nsteps and end_time negative");
	    	
	amrex::IntVect dom_lo(AMREX_D_DECL(0,0,0));
	amrex::IntVect dom_hi(AMREX_D_DECL(ncells[0]-1,ncells[1]-1,ncells[2]-1));
	amrex::Box domain(dom_lo, dom_hi);
	ba.define(domain);
	ba.maxSize(max_grid_size);
	dm.define(ba);
	
	feild_new.define(ba, dm, feild_ncomp,feild_ngrow);
	feild.define(ba, dm, feild_ncomp,feild_ngrow);
	
	amrex::RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},{AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2]}));
	amrex::Array<int,AMREX_SPACEDIM> dom_periodic{AMREX_D_DECL(is_periodic[0],is_periodic[0],is_periodic[0])};
	geom.define(domain,real_box,amrex::CoordSys::cartesian,dom_periodic);
	dx = geom.CellSizeArray();

	set_bc();
	init_prob();
    }


    void set_bc()
    {
	amrex::Vector<amrex::BCRec> bc(feild.nComp());
	for (int n=0; n<feild.nComp(); ++n)
	{
	    for (int idim=0; idim<AMREX_SPACEDIM; ++idim)
	    {
		if (geom.isPeriodic(idim))
		{
		    bc[n].setLo(idim, amrex::BCType::int_dir);
		    bc[n].setHi(idim, amrex::BCType::int_dir);
		}
		else
		{
		    bc[n].setLo(idim, amrex::BCType::foextrap);
		    bc[n].setHi(idim, amrex::BCType::foextrap);
		}
	    }
	}
    }


    void init_prob()
     {
	for (amrex::MFIter mfi(feild); mfi.isValid(); ++mfi)
	{
	    const amrex::Box& box = mfi.validbox();
	    auto const& gnew = feild_new.array(mfi);
	    auto const& g = feild.array(mfi);
	    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
		{
		    g(i,j,k) = 0.5;
		    gnew(i,j,k) = 1;
		});
	}
    }
}
