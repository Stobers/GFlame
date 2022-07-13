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

	amrex::ParmParse pp_prob("prob");
	pp_prob.query("flameloc",flameloc);
	pp_prob.query("flamepert",flamepert);

	
	amrex::IntVect dom_lo(AMREX_D_DECL(0,0,0));
	amrex::IntVect dom_hi(AMREX_D_DECL(ncells[0]-1,ncells[1]-1,ncells[2]-1));
	amrex::Box domain(dom_lo, dom_hi);
	ba.define(domain);
	ba.maxSize(max_grid_size);
	dm.define(ba);
	
	feild_new.define(ba, dm, feild_ncomp,feild_ngrow);
	feild.define(ba, dm, feild_ncomp,feild_ngrow);
	
	amrex::RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},{AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});
	amrex::Array<int,AMREX_SPACEDIM> dom_periodic{AMREX_D_DECL(is_periodic[0],is_periodic[0],is_periodic[0])};
	geom.define(domain,real_box,amrex::CoordSys::cartesian,dom_periodic);
	dx = geom.CellSizeArray();

	set_bc();
	init_gfeild();
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


    void init_gfeild()
    {
	using namespace amrex;
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
	for (MFIter mfi(feild); mfi.isValid(); ++mfi)
	{
	    const Box& box = mfi.validbox();
	    auto const& gnew = feild_new.array(mfi);
	    auto const& g = feild.array(mfi);
	    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
		{
		    AMREX_D_TERM(Real x = prob_lo[0] + i * dx[0];,
				 Real y = prob_lo[1] + j * dx[1];,
				 Real z = prob_lo[2] + k * dx[2]);		 
		    AMREX_D_TERM(Real Lx = prob_hi[0] - prob_lo[0];,
				 Real Ly = prob_hi[1] - prob_lo[1];,
				 Real Lz = prob_hi[2] - prob_lo[2]);
#if (BL_SPACEDIM == 2)
		    Real Ly_perc = Ly / 100;
		    Real loc = flameloc * Ly_perc;
		    Real pert = Ly_perc * (flamepert * 2) *
			(1.0 * std::sin(2 * Pi * 4 * x / Lx) +
			 1.023 * std::sin(2 * Pi * 2 * (x - 0.004598) / Lx) +
			 0.945 * std::sin(2 * Pi * 3 * (x - 0.00712435) / Lx) +
			 1.017 * std::sin(2 * Pi * 5 * (x - 0.0033) / Lx) +
			 0.982 * std::sin(2 * Pi * 5 * (x - 0.014234) / Lx)) / (1 + 1.023 + 0.945 + 1.017 + 0.982);
		    if (y - dx[1] + pert > loc || y + dx[1] + pert < loc)
			g(i,j,k) = (1 / Ly) * (y + pert - loc);
		    else
			g(i,j,k) = 0;
#elif (BL_SPACEDIM == 3)
		    Real Lz_perc = Lz / 100;
		    Real loc = flameloc * Lz_perc;
		    Real pert = Lz_perc * (flamepert * 2) *
			(1.0 * std::sin(2 * Pi * 4 * x / Lx) * std::sin(2 * Pi * 5 * y / Ly) +
			 1.023 * std::sin(2 * Pi * 2 * (x - 0.004598) / Lx) * std::sin(2 * Pi * 4 * (y - 0.0053765) / Ly) +
			 0.945 * std::sin(2 * Pi * 3 * (x - 0.00712435) / Lx) * std::sin(2 * Pi * 3 * (y - 0.02137) / Ly) +
			 1.017 * std::sin(2 * Pi * 5 * (x - 0.0033) / Lx) * std::sin(2 * Pi * 6 * (y - 0.018) / Ly) +
			 0.982 * std::sin(2 * Pi * 5 * (x - 0.014234) / Lx)) / (1 + 1.023 + 0.945 + 1.017 + 0.982);
		    if (z - dx[2] + pert > loc)
			g(i,j,k) = 1;
		    else if (z + dx[2] + pert < loc)
			g(i,j,k) = -1;
		    else
			g(i,j,k) = 0;
#else
		    Abort("invalid ndim");
#endif
		    gnew(i,j,k) = g(i,j,k);
		});
	}
    }
}
