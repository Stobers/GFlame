#include <gflame.hpp>

#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>


namespace Gflame{
    void init()
    {
	amrex::Print() << "initialising domain \n";
	
	amrex::ParmParse pp("dom");
	pp.getarr("ncells",ncells);
	pp.getarr("prob_lo",prob_lo);
	pp.getarr("prob_hi",prob_hi);
	pp.getarr("is_periodic",is_periodic);
	pp.query("max_grid_size",max_grid_size);

	amrex::IntVect dom_lo(AMREX_D_DECL(0,0,0));
	amrex::IntVect dom_hi(AMREX_D_DECL(ncells[0]-1,ncells[1]-1,ncells[2]-1));
	amrex::Box domain(dom_lo, dom_hi);
	ba.define(domain);
	ba.maxSize(max_grid_size);
	dm.define(ba);
	
	int ncomp = 2;
	int ngrow = 0;
	feild.define(ba, dm, ncomp, ngrow);
	
	amrex::RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},{AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2]}));
	amrex::Array<int,AMREX_SPACEDIM> dom_periodic{AMREX_D_DECL(is_periodic[0],is_periodic[0],is_periodic[0])};
	geom.define(domain,real_box,amrex::CoordSys::cartesian,dom_periodic);
	dx = geom.CellSizeArray();

	set_bc();
    }

    
    void set_bc()
    {
	static amrex::Vector<amrex::BCRec> bc(feild.nComp());
	for(int n=0; n<feild.nComp(); ++n)
	{
	    for(int idim=0; idim<AMREX_SPACEDIM; ++idim)
	    {
		if(geom.isPeriodic(idim))
		{
		    amrex::Print() << "is per \n";
		    bc[n].setLo(idim, amrex::BCType::int_dir);
		    bc[n].setHi(idim, amrex::BCType::int_dir);
		}
		else
		{
		    amrex::Print() << "is not per \n";
		    bc[n].setLo(idim, amrex::BCType::foextrap);
		    bc[n].setHi(idim, amrex::BCType::foextrap);
		}
	    }
	}
    }
}
