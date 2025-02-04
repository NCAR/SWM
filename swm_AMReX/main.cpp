/*
 *
 */

#include <numbers>
#include <cmath>

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

#include "swm_mini_app_utils.h"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {

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

    amrex::Print() << "Initial: " << std::endl;
    amrex::Print() << "psi max: " << psi.max(0) << std::endl;
    amrex::Print() << "psi min: " << psi.min(0) << std::endl;
    amrex::Print() << "p max: " << p.max(0) << std::endl;
    amrex::Print() << "p min: " << p.min(0) << std::endl;
    amrex::Print() << "u max: " << u.max(0) << std::endl;
    amrex::Print() << "u min: " << u.min(0) << std::endl;
    amrex::Print() << "v max: " << v.max(0) << std::endl;
    amrex::Print() << "v min: " << v.min(0) << std::endl;

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

    // Constants used in time stepping loop
    const double fsdx = 4.0/dx;
    const double fsdy = 4.0/dy;
    double tdt = dt;
    const double alpha = 0.001; 

    for (int time_step = 0; time_step < n_time_steps; ++time_step)
    {
        // Sets: cu, cv, h, z
        UpdateIntermediateVariables(fsdx, fsdy, geom,
                                     p, u, v,
                                     cu, cv, h, z);


        // Sets: p_new, u_new, v_new
        UpdateNewVariables(dx, dy, tdt, geom,
                           p_old, u_old, v_old, cu, cv, h, z,
                           p_new, u_new, v_new);


        // Sets: p_old, u_old, v_old
        UpdateOldVariables(alpha, time_step, geom,
                           p, u, v,
                           p_new, u_new, v_new,
                           p_old, u_old, v_old);

        // Sets: p, u, v
        UpdateVariables(geom, u_new, v_new, p_new, u, v, p);

        if (time_step == 0) {
            tdt = tdt + tdt;
        }

        time = time + dt;

        // Write a plotfile of the current data (plot_interval was defined in the inputs file)
        if (plot_interval > 0 && time_step%plot_interval == 0)
        {
            WriteOutput(psi, p, u, v, geom, time, time_step, output_values);
        }

    }

    amrex::Print() << "Final: " << std::endl;
    amrex::Print() << "p max: " << p.max(0) << std::endl;
    amrex::Print() << "p min: " << p.min(0) << std::endl;
    amrex::Print() << "u max: " << u.max(0) << std::endl;
    amrex::Print() << "u min: " << u.min(0) << std::endl;
    amrex::Print() << "v max: " << v.max(0) << std::endl;
    amrex::Print() << "v min: " << v.min(0) << std::endl;

    }
    amrex::Finalize();
    return 0;
}