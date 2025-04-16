#include <string>
#include <cmath>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

#include "swm_mini_app_utils.h"
#include "swm_mini_app_kernels.h"

void ParseInput(int & nx, int & ny,
                amrex::Real & dx, amrex::Real & dy,
                int & max_chunk_size,
                int & n_time_steps, amrex::Real & dt,
                int & plot_interval)
{
    // ParmParse is way of reading inputs from the inputs file
    // pp.get means we require the inputs file to have it
    // pp.query means we optionally need the inputs file to have it - but we must supply a default here
    amrex::ParmParse pp;

    pp.get("nx",nx);
    pp.get("ny",ny);

    pp.get("dx",dx);
    pp.get("dy",dy);

    pp.get("max_chunk_size",max_chunk_size);

    pp.get("n_time_steps",n_time_steps);

    pp.get("dt",dt);

    // Default plot_interval to -1, allow us to set it to something else in the inputs file
    //  If plot_interval < 0 then no plot files will be written
    plot_interval = -1;
    pp.query("plot_interval",plot_interval);

    return;
}

void DefineCellCenteredMultiFab(const int nx, const int ny,
                                const int max_chunk_size,
                                amrex::MultiFab & cell_centered_MultiFab)
{
    // lower and upper indices of domain
    const amrex::IntVect domain_low_index(0,0);
    const amrex::IntVect domain_high_index(nx-1, ny-1);
    
    // create box of indicies for cells
    const amrex::Box cell_centered_box(domain_low_index, domain_high_index);

    // initialize the boxarray "cell_box_array" from the single box "cell_centered_box"
    amrex::BoxArray cell_box_array(cell_centered_box);
    //cell_box_array.define(cell_centered_box);

    // break up boxarray "cell_box_array" into chunks no larger than "max_chunk_size" along a direction
    cell_box_array.maxSize(max_chunk_size);

    // assigns processor to each box in the box array
    amrex::DistributionMapping distribution_mapping(cell_box_array);

    //amrex::Print() << "max_chunk_size: " << max_chunk_size << std::endl;
    //amrex::Print() << "cell_box_array: " << cell_box_array << std::endl;
    //amrex::Print() << "distribution mapping: " << distribution_mapping << std::endl;

    // number of components for each array
    int Ncomp = 1;

    // number of ghost cells for each array
    int Nghost = 1;

    cell_centered_MultiFab.define(cell_box_array, distribution_mapping, Ncomp, Nghost);

    return;
}

void DefineXFaceMultiFab(const amrex::MultiFab & cell_centered_MultiFab,
                         amrex::MultiFab & x_face_MultiFab)
{
    AMREX_ASSERT(cell_centered_MultiFab.is_cell_centered());

    const amrex::BoxArray cell_box_array = cell_centered_MultiFab.boxArray();

    const amrex::BoxArray x_face_box_array = amrex::convert(cell_box_array, {1,0});

    x_face_MultiFab.define(x_face_box_array, cell_centered_MultiFab.DistributionMap(), cell_centered_MultiFab.nComp(), cell_centered_MultiFab.nGrow());

    return;
}

void DefineYFaceMultiFab(const amrex::MultiFab & cell_centered_MultiFab,
                         amrex::MultiFab & y_face_MultiFab)
{
    AMREX_ASSERT(cell_centered_MultiFab.is_cell_centered());

    const amrex::BoxArray cell_box_array = cell_centered_MultiFab.boxArray();

    const amrex::BoxArray y_face_box_array = amrex::convert(cell_box_array, {0,1});

    y_face_MultiFab.define(y_face_box_array, cell_centered_MultiFab.DistributionMap(), cell_centered_MultiFab.nComp(), cell_centered_MultiFab.nGrow());

    return;
}

void DefineNodalMultiFab(const amrex::MultiFab & cell_centered_MultiFab,
                         amrex::MultiFab & nodal_MultiFab)
{
    AMREX_ASSERT(cell_centered_MultiFab.is_cell_centered());

    const amrex::BoxArray cell_box_array = cell_centered_MultiFab.boxArray();


    amrex::BoxArray surrounding_nodes_box_array = cell_box_array;
    surrounding_nodes_box_array.surroundingNodes();

    nodal_MultiFab.define(surrounding_nodes_box_array, cell_centered_MultiFab.DistributionMap(), cell_centered_MultiFab.nComp(), cell_centered_MultiFab.nGrow());

    return;
}

void InitializeGeometry(const int nx, const int ny,
                        const amrex::Real dx, const amrex::Real dy,
                        amrex::Geometry & geom)
{
  // lower and upper indices of domain
  const amrex::IntVect domain_low_index(0,0);
  const amrex::IntVect domain_high_index(nx-1, ny-1);

  // create box of indicies for cells
  const amrex::Box cell_centered_box(domain_low_index, domain_high_index);

  // physical min and max boundaries of cells
  const amrex::RealBox real_box({0, 0},
                                {nx*dx, ny*dy});

  // This, a value of 0, says we are using Cartesian coordinates
  int coordinate_system = 0;

  // This sets the boundary conditions in each direction to periodic
  amrex::Array<int,AMREX_SPACEDIM> is_periodic {1,1};

  // This defines a Geometry object
  geom.define(cell_centered_box, real_box, coordinate_system, is_periodic);
 // geom.define(cell_centered_box, real_box, amrex::CoordSys::cartesian, is_periodic); // Could use an amrex defined enum instead of an int to specify the coordinate system

  return;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real LinearMapCoordinates(const amrex::Real x, 
                                 const amrex::Real x_min, const amrex::Real x_max,
                                 const amrex::Real xi_min, const amrex::Real xi_max)
{
    return x_min + ((xi_max-xi_min)/(x_max-x_min))*x;
}

void InitializeVariables(const amrex::Geometry & geom,
                         amrex::MultiFab & psi,
                         amrex::MultiFab & p,
                         amrex::MultiFab & u,
                         amrex::MultiFab & v)
{

    const amrex::Real x_min = geom.ProbLo(0);
    const amrex::Real x_max = geom.ProbHi(0);
    const amrex::Real y_min = geom.ProbLo(1);
    const amrex::Real y_max = geom.ProbHi(1);

    const amrex::Real dx = geom.CellSize(0);
    const amrex::Real dy = geom.CellSize(1);

    ////////////////////////////////////////////////////////////////////////// 
    // Initialization of stream function (psi)
    ////////////////////////////////////////////////////////////////////////// 

    // coefficient for initialization psi
    const amrex::Real a = 1000000;
    const double pi = 4. * std::atan(1.);

    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& phi_array = psi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const amrex::Real x_cell_center = (i+0.5) * dx;
            const amrex::Real y_cell_center = (j+0.5) * dy;

            const amrex::Real x_transformed = LinearMapCoordinates(x_cell_center, x_min, x_max, 0.0, 2*pi);
            const amrex::Real y_transformed = LinearMapCoordinates(y_cell_center, y_min, y_max, 0.0, 2*pi);

            phi_array(i,j,k) = a*std::sin(x_transformed)*std::sin(y_transformed);
        });
    }
    
    psi.FillBoundary(geom.periodicity());

    ////////////////////////////////////////////////////////////////////////// 
    // Initialization of pressure (p)
    ////////////////////////////////////////////////////////////////////////// 

    // coefficient for pressure
    double el = geom.ProbLength(0);
    amrex::Real pcf = (pi * pi * a * a)/(el * el);

    for (amrex::MFIter mfi(p); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& p_array = p.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            amrex::Real x_node = i * dx;
            amrex::Real y_node = j * dy;
            
            const amrex::Real x_transformed = LinearMapCoordinates(x_node, x_min, x_max, 0.0, 2*pi);
            const amrex::Real y_transformed = LinearMapCoordinates(y_node, y_min, y_max, 0.0, 2*pi);

            p_array(i,j,k) = pcf * (std::cos(2*x_transformed) + std::cos(2*y_transformed)) + 5000;
        });
    }

    p.FillBoundary(geom.periodicity());

    ////////////////////////////////////////////////////////////////////////// 
    // Initialization of x velocity (u)
    ////////////////////////////////////////////////////////////////////////// 

    for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& u_array = u.array(mfi);
        const amrex::Array4<amrex::Real const>& phi_old_array = psi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            u_array(i,j,k) = -(phi_old_array(i,j,k)-phi_old_array(i,j-1,k))/dy;
        });
    }

    u.FillBoundary(geom.periodicity());

    ////////////////////////////////////////////////////////////////////////// 
    // Initialization of y velocity (v)
    ////////////////////////////////////////////////////////////////////////// 

    for (amrex::MFIter mfi(v); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& v_array = v.array(mfi);
        const amrex::Array4<amrex::Real const>& phi_old_array = psi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            v_array(i,j,k) = (phi_old_array(i,j,k)-phi_old_array(i-1,j,k))/dx;
        });
    }

    v.FillBoundary(geom.periodicity());

    return;
}


void WriteOutput(const amrex::MultiFab & psi,
                 const amrex::MultiFab & p,
                 const amrex::MultiFab & u,
                 const amrex::MultiFab & v,
                 const amrex::Geometry & geom,
                 const amrex::Real time,
                 const int time_step,
                 amrex::MultiFab & output_values)
{

    // Interpolate all values to cell centers
    for (amrex::MFIter mfi(output_values); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real const>& phi_old_array = psi.const_array(mfi);
        const amrex::Array4<amrex::Real const>& p_array = p.const_array(mfi);
        const amrex::Array4<amrex::Real const>& u_array = u.const_array(mfi);
        const amrex::Array4<amrex::Real const>& v_array = v.const_array(mfi);

        const amrex::Array4<amrex::Real>& output_values_array = output_values.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            output_values_array(i,j,k,0) = phi_old_array(i,j,k);
            output_values_array(i,j,k,1) = (p_array(i,j,k) + p_array(i+1,j,k) + p_array(i,j+1,k) + p_array(i+1,j+1,k))/4.0;
            output_values_array(i,j,k,2) = (u_array(i,j,k) + u_array(i,j+1,k))/2.0;
            output_values_array(i,j,k,3) = (v_array(i,j,k) + v_array(i+1,j,k))/2.0;
        });
    }

    // Write output file
    int min_digits = 5;
    const std::string& pltfile = amrex::Concatenate("plt", time_step, min_digits);
    amrex::WriteSingleLevelPlotfile(pltfile, output_values, {"psi", "p", "u", "v"}, geom, time, time_step);

    return;
}

amrex::MultiFab CreateMultiFab(const amrex::MultiFab & mf)
{
    return amrex::MultiFab(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrow());
}

void Copy(const amrex::MultiFab & src, amrex::MultiFab & dest)
{
    // TODO: Error check that the number of components and ghost cells are the same. Maybe also all the other properties of the multifab

    const int src_starting_component_index = 0;
    const int dest_starting_component_index = 0;
    amrex::MultiFab::Copy(dest, src, src_starting_component_index, dest_starting_component_index, src.nComp(), src.nGrow());
    return;
}

void Swap(amrex::MultiFab & src, amrex::MultiFab & dest)
{
    // TODO: Error check that the number of components and ghost cells are the same. Maybe also all the other properties of the multifab 

    const int src_starting_component_index = 0;
    const int dest_starting_component_index = 0;
    amrex::MultiFab::Swap(dest, src, src_starting_component_index, dest_starting_component_index, src.nComp(), src.nGrow());
    return;
}

void UpdateIntermediateVariables(amrex::Real dx, amrex::Real dy, const amrex::Geometry& geom,
                                 const amrex::MultiFab& p, const amrex::MultiFab& u, const amrex::MultiFab& v,
                                 amrex::MultiFab& cu, amrex::MultiFab& cv, amrex::MultiFab& h, amrex::MultiFab& z)
{
    BL_PROFILE("UpdateIntermediateVariables()");

    const double fsdx = 4.0/dx;
    const double fsdy = 4.0/dy;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(p, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();

        // Read only arrays
        const amrex::Array4<amrex::Real const>& p_array = p.const_array(mfi);
        const amrex::Array4<amrex::Real const>& u_array = u.const_array(mfi);
        const amrex::Array4<amrex::Real const>& v_array = v.const_array(mfi);

        // Write arrays
        const amrex::Array4<amrex::Real>& cu_array = cu.array(mfi);
        const amrex::Array4<amrex::Real>& cv_array = cv.array(mfi);
        const amrex::Array4<amrex::Real>& h_array = h.array(mfi);
        const amrex::Array4<amrex::Real>& z_array = z.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            UpdateIntermediateVariablesKernel(i, j, k, fsdx, fsdy,
                                              p_array, u_array, v_array,
                                              cu_array, cv_array, h_array, z_array);
        });
    }

    cu.FillBoundary(geom.periodicity());
    cv.FillBoundary(geom.periodicity());
    h.FillBoundary(geom.periodicity());
    z.FillBoundary(geom.periodicity());

    return;
}

void UpdateNewVariables(const double dx, const double dy, const double tdt, const amrex::Geometry& geom,
                        const amrex::MultiFab& p_old, const amrex::MultiFab& u_old, const amrex::MultiFab& v_old,
                        const amrex::MultiFab& cu, const amrex::MultiFab& cv, const amrex::MultiFab& h, const amrex::MultiFab& z,
                        amrex::MultiFab& p_new, amrex::MultiFab& u_new, amrex::MultiFab& v_new)
{
    BL_PROFILE("UpdateNewVariables()");

    // defined here because tdt changes after first time step
    const double tdtsdx = tdt / dx;
    const double tdtsdy = tdt / dy;
    const double tdts8  = tdt / 8.0;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(p_old, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();

        // Read only arrays
        const amrex::Array4<amrex::Real const>& p_old_array = p_old.const_array(mfi);
        const amrex::Array4<amrex::Real const>& u_old_array = u_old.const_array(mfi);
        const amrex::Array4<amrex::Real const>& v_old_array = v_old.const_array(mfi);
        const amrex::Array4<amrex::Real const>& cu_array = cu.const_array(mfi);
        const amrex::Array4<amrex::Real const>& cv_array = cv.const_array(mfi);
        const amrex::Array4<amrex::Real const>& h_array = h.const_array(mfi);
        const amrex::Array4<amrex::Real const>& z_array = z.const_array(mfi);

        // Write arrays
        const amrex::Array4<amrex::Real>& p_new_array = p_new.array(mfi);
        const amrex::Array4<amrex::Real>& u_new_array = u_new.array(mfi);
        const amrex::Array4<amrex::Real>& v_new_array = v_new.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            UpdateNewVariablesKernel(i, j, k, 
                                     tdtsdx, tdtsdy, tdts8,
                                     p_old_array, u_old_array, v_old_array,
                                     cu_array, cv_array, h_array, z_array,
                                     p_new_array, u_new_array, v_new_array);

        });
    }

    return;
}

void UpdateOldVariables(const double alpha, const int time_step, const amrex::Geometry& geom,
                        const amrex::MultiFab& p, const amrex::MultiFab& u, const amrex::MultiFab& v, 
                        const amrex::MultiFab& p_new, const amrex::MultiFab& u_new, const amrex::MultiFab& v_new, 
                        amrex::MultiFab& p_old, amrex::MultiFab& u_old, amrex::MultiFab& v_old)
{
    BL_PROFILE("UpdateOldVariables()");
    if (time_step > 0) {

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(p, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();

            // Read only arrays
            const amrex::Array4<amrex::Real const>& p_array = p.const_array(mfi);
            const amrex::Array4<amrex::Real const>& u_array = u.const_array(mfi);
            const amrex::Array4<amrex::Real const>& v_array = v.const_array(mfi);
            const amrex::Array4<amrex::Real const>& p_new_array = p_new.const_array(mfi);
            const amrex::Array4<amrex::Real const>& u_new_array = u_new.const_array(mfi);
            const amrex::Array4<amrex::Real const>& v_new_array = v_new.const_array(mfi);

            // Write arrays
            const amrex::Array4<amrex::Real>& p_old_array = p_old.array(mfi);
            const amrex::Array4<amrex::Real>& u_old_array = u_old.array(mfi);
            const amrex::Array4<amrex::Real>& v_old_array = v_old.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                UpdateOldVariablesKernel(i, j, k, 
                                         alpha,
                                         p_array, u_array, v_array,
                                         p_new_array, u_new_array, v_new_array,
                                         p_old_array, u_old_array, v_old_array);
            });
        }
    } else {

        Copy(u, u_old);
        Copy(v, v_old);
        Copy(p, p_old);
    }

    return;
}

void UpdateVariables(const amrex::Geometry& geom, 
                     amrex::MultiFab& u_new, amrex::MultiFab& v_new, amrex::MultiFab& p_new,
                     amrex::MultiFab& u, amrex::MultiFab& v, amrex::MultiFab& p)
{
    BL_PROFILE("UpdateVariables()");

    Swap(u_new, u);
    Swap(v_new, v);
    Swap(p_new, p);

    u.FillBoundary(geom.periodicity());
    v.FillBoundary(geom.periodicity());
    p.FillBoundary(geom.periodicity());

    return;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateXFaceToNode(const amrex::Array4<amrex::Real const>& phi_x_face, 
                                   const int i, const int j, const int k)
{
    amrex::Real phi_node = 0.5*(phi_x_face(i,j-1,k) + phi_x_face(i,j,k));
    return phi_node;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateXFaceToYFace(const amrex::Array4<amrex::Real const>& phi_x_face, 
                                    const int i, const int j, const int k)
{
    amrex::Real phi_y_face = 0.25*(phi_x_face(i,j-1,k) + phi_x_face(i+1,j-1,k) + phi_x_face(i,j,k) + phi_x_face(i+1,j,k) );
    return phi_y_face;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateXFaceToCellCenter(const amrex::Array4<amrex::Real const>& phi_x_face, 
                                         const int i, const int j, const int k)
{
    amrex::Real phi_cell_center = 0.5*(phi_x_face(i,j,k) + phi_x_face(i+1,j,k) );
    return phi_cell_center;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateYFaceToNode(const amrex::Array4<amrex::Real const>& phi_y_face, 
                                   const int i, const int j, const int k)
{
    amrex::Real phi_node = 0.5*(phi_y_face(i-1,j,k) + phi_y_face(i,j,k));
    return phi_node;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateYFaceToXFace(const amrex::Array4<amrex::Real const>& phi_y_face, 
                                    const int i, const int j, const int k)
{
    amrex::Real phi_x_face = 0.25*(phi_y_face(i-1,j,k) + phi_y_face(i,j,k) + phi_y_face(i-1,j+1,k) + phi_y_face(i,j+1,k) );
    return phi_x_face;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateYFaceToCellCenter(const amrex::Array4<amrex::Real const>& phi_y_face, 
                                         const int i, const int j, const int k)
{
    amrex::Real phi_cell_center = 0.5*(phi_y_face(i,j,k) + phi_y_face(i,j+1,k));
    return phi_cell_center;
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateNodeToXFace(const amrex::Array4<amrex::Real const>& phi_node, 
                                   const int i, const int j, const int k)
{
    amrex::Real phi_x_face = 0.5*(phi_node(i,j,k) + phi_node(i,j+1,k));
    return phi_x_face;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateNodeToYFace(const amrex::Array4<amrex::Real const>& phi_node, 
                                   const int i, const int j, const int k)
{
    amrex::Real phi_y_face = 0.5*(phi_node(i,j,k) + phi_node(i+1,j,k));
    return phi_y_face;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateNodeToCellCenter(const amrex::Array4<amrex::Real const>& phi_node, 
                                        const int i, const int j, const int k)
{
    amrex::Real phi_cell_center = 0.25*(phi_node(i,j,k) + phi_node(i+1,j,k) + phi_node(i,j+1,k) + phi_node(i+1,j+1,k));
    return phi_cell_center;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateCellCenterToNode(const amrex::Array4<amrex::Real const>& phi_cell_center, 
                                        const int i, const int j, const int k)
{
    amrex::Real phi_node = 0.25*(phi_cell_center(i-1,j-1,k) + phi_cell_center(i,j-1,k) + phi_cell_center(i-1,j,k) + phi_cell_center(i,j,k));
    return phi_node;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateCellCenterToXFace(const amrex::Array4<amrex::Real const>& phi_cell_center, 
                                        const int i, const int j, const int k)
{
    amrex::Real phi_x_face = 0.5*(phi_cell_center(i-1,j,k) + phi_cell_center(i,j,k));
    return phi_x_face;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateCellCenterToYFace(const amrex::Array4<amrex::Real const>& phi_cell_center, 
                                        const int i, const int j, const int k)
{
    amrex::Real phi_y_face = 0.5*(phi_cell_center(i,j-1,k) + phi_cell_center(i,j,k));
    return phi_y_face;
}


void UpdateNewVariables(const double dx, const double dy, const double dt,
                        const amrex::MultiFab& p_old, const amrex::MultiFab& u_old, const amrex::MultiFab& v_old,
                        amrex::MultiFab& p_new, amrex::MultiFab& u_new, amrex::MultiFab& v_new)
{
    BL_PROFILE("UpdateNewVariables()");

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(p_old, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();

        // Read only arrays
        const amrex::Array4<amrex::Real const>& p_old_array = p_old.const_array(mfi);
        const amrex::Array4<amrex::Real const>& u_old_array = u_old.const_array(mfi);
        const amrex::Array4<amrex::Real const>& v_old_array = v_old.const_array(mfi);

        // Write arrays
        const amrex::Array4<amrex::Real>& p_new_array = p_new.array(mfi);
        const amrex::Array4<amrex::Real>& u_new_array = u_new.array(mfi);
        const amrex::Array4<amrex::Real>& v_new_array = v_new.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            UpdateNewVariablesKernel(i, j, k, 
                                     dx, dy, dt,
                                     p_old_array, u_old_array, v_old_array,
                                     p_new_array, u_new_array, v_new_array);

        });
    }

    return;
}