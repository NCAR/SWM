#include <string>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

#include "swm_mini_app_utils.h"
#include "swm_mini_app_kernels.h"

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