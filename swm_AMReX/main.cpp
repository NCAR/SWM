/*
 *
 */

#include <cmath>

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

#include "swm_mini_app_utils.h"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
    BL_PROFILE("main()");

    // ***********************************************************************
    // Simulation Parameters Set Via Input File
    // ***********************************************************************

    // Number of cells on each direction
    int nx;
    int ny;

    // Cell size in each direction
    amrex::Real dx;
    amrex::Real dy;

    // Mesh will be broken into chunks of up to max_chunk_size
    int max_chunk_size;

    // Number of time steps to take
    int n_time_steps;

    // Size of time step
    amrex::Real dt;

    // How often to write a plotfile
    //     Optional argument. If left out no plot files will be written.
    int plot_interval;

    // Set parameter values from inputs file
    ParseInput(nx, ny, dx, dy, max_chunk_size,
               n_time_steps, dt, plot_interval);

    // ***********************************************************************
    // Define arrays
    // ***********************************************************************

    amrex::MultiFab psi;
    DefineCellCenteredMultiFab(nx, ny, max_chunk_size, psi);

    amrex::MultiFab v;
    DefineXFaceMultiFab(psi, v);

    amrex::MultiFab u;
    DefineYFaceMultiFab(psi, u);

    amrex::MultiFab p;
    DefineNodalMultiFab(psi, p);
    
    // **********************************
    // Initialize Data
    // **********************************

    // AMReX object to hold domain meta data... Like the physical size of the domain and if it is periodic in each direction
    amrex::Geometry geom;
    InitializeGeometry(nx, ny, dx, dy, geom);

    InitializeVariables(geom, psi, p, u, v);

    //amrex::Print() << "Initial: " << std::endl;
    //amrex::Print() << "psi max: " << psi.max(0) << std::endl;
    //amrex::Print() << "psi min: " << psi.min(0) << std::endl;
    //amrex::Print() << "p max: " << p.max(0) << std::endl;
    //amrex::Print() << "p min: " << p.min(0) << std::endl;
    //amrex::Print() << "u max: " << u.max(0) << std::endl;
    //amrex::Print() << "u min: " << u.min(0) << std::endl;
    //amrex::Print() << "v max: " << v.max(0) << std::endl;
    //amrex::Print() << "v min: " << v.min(0) << std::endl;

    // **********************************
    // Write initial plot file
    // **********************************

    amrex::Real time = 0.0;

    // Interpolate the values to the cell center for writing output
    amrex::MultiFab output_values(psi.boxArray(), psi.DistributionMap(), 4, 0);

    if (plot_interval > 0)
    {
        int time_step = 0;
        WriteOutput(psi, p, u, v, geom, time, time_step, output_values);
    }

    // **********************************************
    // Intermediate Values used in time stepping loop
    // **********************************************

    // The product of pressure and x-velocity. 
    // Called cu for consistency with the other version of the mini-app
    // Stored on the y faces (same locations as u)
    amrex::MultiFab cu = CreateMultiFab(u);

    // The product of pressure and y-velocity. 
    // Called cv for consistency with the other version of the mini-app
    // Stored on the x faces (same locations as v)
    amrex::MultiFab cv = CreateMultiFab(v);

    // The potential vorticity. 
    // Called z for consistency with the other version of the mini-app
    // Stored on the cell centers (same locations as psi)
    amrex::MultiFab z = CreateMultiFab(psi);

    // The term (P + 1/2(V dot V)). The gradient of this term appears on the right hand side of the momentum equations.
    // Called h for consistency with the other version of the mini-app
    // Stored on the nodal points (same locations as p)
    amrex::MultiFab h = CreateMultiFab(p);

    // Arrays to hold the primary variables (u,v,p) that are updated when time stepping
    amrex::MultiFab u_old = CreateMultiFab(u);
    amrex::MultiFab v_old = CreateMultiFab(v);
    amrex::MultiFab p_old = CreateMultiFab(p);
    amrex::MultiFab u_new = CreateMultiFab(u);
    amrex::MultiFab v_new = CreateMultiFab(v);
    amrex::MultiFab p_new = CreateMultiFab(p);

    // For the first time step the {u,v,p}_old values are initialized to match {u,v,p}.
    Copy(u, u_old);
    Copy(v, v_old);
    Copy(p, p_old);


    auto printMultiFabInfo = [](const amrex::MultiFab& multiFab_to_print) {
        for (int i = 0; i < multiFab_to_print.boxArray().size(); ++i) {
            amrex::Print() << "box " << std::to_string(i) << " size " << multiFab_to_print.boxArray()[i].size() << std::endl;
        }
        amrex::Print() << "distribution map " << multiFab_to_print.DistributionMap() << std::endl;
    };

    // Store the MultiFab objects and their names in a vector of pairs. This will allow us to iterate over them and call the lambda function.
    std::vector<std::pair<std::string, const amrex::MultiFab&>> multiFabs = {
        {"psi", psi},
        {"u", u},
        {"v", v},
        {"p", p},
        {"cu", cu},
        {"cv", cv},
        {"h", h},
        {"z", z},
        {"u_new", u_new},
        {"v_new", v_new},
        {"p_new", p_new},
        {"u_old", u_old},
        {"v_old", v_old},
        {"p_old", p_old}
    };

    // Iterate over the arrays and call the lambda function
    for (const auto& [name, multiFab] : multiFabs) {
        amrex::Print() << name << ": " << std::endl;
        printMultiFabInfo(multiFab);
    }

    // Constants used in time stepping loop
    double tdt = dt;
    const double alpha = 0.001; 
    BL_PROFILE_VAR_NS("SWM:Loop 100",loop100);
    BL_PROFILE_VAR_NS("SWM:Loop 200",loop200);
    BL_PROFILE_VAR_NS("SWM:Loop 300",loop300);
    BL_PROFILE_VAR_NS("SWM:Total",total);

    BL_PROFILE_VAR_START(total);
    for (int time_step = 1; time_step <= n_time_steps; ++time_step)
    {
	BL_PROFILE_VAR_START(loop100);
        // Sets: cu, cv, h, z
        UpdateIntermediateVariables(dx, dy, geom,
                                     p, u, v,
                                     cu, cv, h, z);
	BL_PROFILE_VAR_STOP(loop100);


        // Sets: p_new, u_new, v_new
	BL_PROFILE_VAR_START(loop200);
        UpdateNewVariables(dx, dy, tdt, geom,
                           p_old, u_old, v_old, cu, cv, h, z,
                           p_new, u_new, v_new);
	BL_PROFILE_VAR_STOP(loop200);


	BL_PROFILE_VAR_START(loop300);
        // Sets: p_old, u_old, v_old
        UpdateOldVariables(alpha, time_step, geom,
                           p, u, v,
                           p_new, u_new, v_new,
                           p_old, u_old, v_old);

        // Sets: p, u, v
        UpdateVariables(geom, u_new, v_new, p_new, u, v, p);
	BL_PROFILE_VAR_STOP(loop300);

        time = time + dt;

        // Write a plotfile of the current data (plot_interval was defined in the inputs file)
        if (plot_interval > 0 && time_step%plot_interval == 0)
        {
            WriteOutput(psi, p, u, v, geom, time, time_step, output_values);
        }

        if (time_step == 0) {
            tdt = tdt + tdt;
        }

    }

    //amrex::Print() << "Final: " << std::endl;
    //amrex::Print() << "p max: " << p.max(0) << std::endl;
    //amrex::Print() << "p min: " << p.min(0) << std::endl;
    //amrex::Print() << "u max: " << u.max(0) << std::endl;
    //amrex::Print() << "u min: " << u.min(0) << std::endl;
    //amrex::Print() << "v max: " << v.max(0) << std::endl;
    //amrex::Print() << "v min: " << v.min(0) << std::endl;
    BL_PROFILE_VAR_STOP(total);

    }

    amrex::Finalize();
    return 0;
}
