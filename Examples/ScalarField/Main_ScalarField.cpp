/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// AMReX includes
#include <AMReX.H> //Gives us amrex::Print()
#include <AMReX_ParmParse.H>

// Our general includes
#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "ScalarFieldLevel.hpp"

// Chombo namespace
using namespace amrex;

int runGRTeclyn(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.

    GRParmParse pp;
    SimulationParameters sim_params(pp);

    if (sim_params.just_check_params)
    {
        return 0;
    }

    // The line below selects the problem that is simulated
    // (To simulate a different problem, define a new child of AMRLevel
    // and an associated LevelFactory)

    GRAMR::set_simulation_parameters(sim_params);

    DefaultLevelFactory<ScalarFieldLevel> scalar_field_level_bld;

    GRAMR gr_amr(&scalar_field_level_bld); // check that this isn't meant to be
                                           // a special class for scalar fields
    gr_amr.init(0., sim_params.stop_time);

    while (
        (gr_amr.okToContinue() != 0) &&
        (gr_amr.levelSteps(0) < sim_params.max_steps ||
         sim_params.max_steps < 0) &&
        (gr_amr.cumTime() < sim_params.stop_time || sim_params.stop_time < 0.0))
    {
        gr_amr.coarseTimeStep(sim_params.stop_time);
    }

    // Write final checkpoint and plotfile
    if (gr_amr.stepOfLastCheckPoint() < gr_amr.levelSteps(0) &&
        sim_params.checkpoint_interval >= 0)
    {
        gr_amr.checkPoint();
    }

    if (gr_amr.stepOfLastPlotFile() < gr_amr.levelSteps(0) &&
        sim_params.plot_interval >= 0)
    {
        gr_amr.writePlotFile();
    }

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRTeclyn(argc, argv);

    if (status == 0)
        amrex::Print() << "GRChombo finished." << std::endl;
    else
        amrex::Print() << "GRChombo failed with return code " << status
                       << std::endl;

    mainFinalize();
    return status;
}
