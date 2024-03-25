#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include "KleinGordon.hpp"

using namespace amrex;

amrex::LevelBld* getLevelBld ();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    int  max_step = -1;
    Real strt_time = 0.0;
    Real stop_time = -1.0;

   
    {
        ParmParse pp;
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);


    }
    if (max_step < 0 && stop_time < 0.0) {
        amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    {
        auto amr = std::make_unique<KleinGordon>(getLevelBld());

        amr->init(strt_time,stop_time);

	
        while ( amr->okToContinue() &&
                (amr->levelSteps(0) < max_step || max_step < 0) &&
                (amr->cumTime() < stop_time || stop_time < 0.0) )
        {
            amr->coarseTimeStep(stop_time);
	    int current_step = amr->levelSteps(0);

	    
	    //I think this gives the current step for level 0
	    if (amr->plotInt() > 0 && current_step%amr->plotInt()==0)
	      {
		const std::string& output_fname = amrex::Concatenate("plt", current_step, 5);
		//		amrex::Print() << AMREX_SPACEDIM << "\n\n";
		//		WriteMultiLevelPlotfile(output_fname, amr.max_levels, 
	      }
	    
        }

        // Write final checkpoint and plotfile

        if (amr->stepOfLastCheckPoint() < amr->levelSteps(0)) {
            amr->checkPoint();
        }
        if (amr->stepOfLastPlotFile() < amr->levelSteps(0)) {
            amr->writePlotFile();
        }
    }

    amrex::Finalize();
}
