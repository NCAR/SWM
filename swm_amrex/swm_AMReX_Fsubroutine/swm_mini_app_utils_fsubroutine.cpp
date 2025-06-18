#include <string>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

#include "swm_mini_app_utils.h"
#include "funcF.H"

void UpdateIntermediateVariables(amrex::Real  dx, amrex::Real  dy, const amrex::Geometry& geom,
                                 const amrex::MultiFab& p, const amrex::MultiFab& u, const amrex::MultiFab& v,
                                 amrex::MultiFab& cu, amrex::MultiFab& cv, amrex::MultiFab& h, amrex::MultiFab& z)
{

    amrex::Real fsdx = 4.0/dx;
    amrex::Real fsdy = 4.0/dy;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(p, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();

	UpdateIntermediateVariablesSub( BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_ANYD(p[mfi]),
                     BL_TO_FORTRAN_ANYD(u[mfi]),
                     BL_TO_FORTRAN_ANYD(v[mfi]),
                     BL_TO_FORTRAN_ANYD(cu[mfi]),
                     BL_TO_FORTRAN_ANYD(cv[mfi]),
                     BL_TO_FORTRAN_ANYD(h[mfi]),
                     BL_TO_FORTRAN_ANYD(z[mfi]),
                     fsdx,fsdy);
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

	UpdateNewVariablesSub( BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_ANYD(p_old[mfi]),
                     BL_TO_FORTRAN_ANYD(u_old[mfi]),
                     BL_TO_FORTRAN_ANYD(v_old[mfi]),
                     BL_TO_FORTRAN_ANYD(cu[mfi]),
                     BL_TO_FORTRAN_ANYD(cv[mfi]),
                     BL_TO_FORTRAN_ANYD(h[mfi]),
                     BL_TO_FORTRAN_ANYD(z[mfi]),
                     BL_TO_FORTRAN_ANYD(p_new[mfi]),
                     BL_TO_FORTRAN_ANYD(u_new[mfi]),
                     BL_TO_FORTRAN_ANYD(v_new[mfi]),
                     tdtsdx,tdtsdy,tdts8);
    }

    return;
}

void UpdateOldVariables(const double alpha, const int time_step, const amrex::Geometry& geom,
                        const amrex::MultiFab& p, const amrex::MultiFab& u, const amrex::MultiFab& v, 
                        const amrex::MultiFab& p_new, const amrex::MultiFab& u_new, const amrex::MultiFab& v_new, 
                        amrex::MultiFab& p_old, amrex::MultiFab& u_old, amrex::MultiFab& v_old)
{
    if (time_step > 0) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
	for (amrex::MFIter mfi(p, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
	    UpdateOldVariablesSub(BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_ANYD(p_new[mfi]),
                     BL_TO_FORTRAN_ANYD(u_new[mfi]),
                     BL_TO_FORTRAN_ANYD(v_new[mfi]),
                     BL_TO_FORTRAN_ANYD(p[mfi]),
                     BL_TO_FORTRAN_ANYD(u[mfi]),
                     BL_TO_FORTRAN_ANYD(v[mfi]),
                     BL_TO_FORTRAN_ANYD(p_old[mfi]),
                     BL_TO_FORTRAN_ANYD(u_old[mfi]),
                     BL_TO_FORTRAN_ANYD(v_old[mfi]),
                     alpha);
        }
    } else {

        Copy(u, u_old);
        Copy(v, v_old);
        Copy(p, p_old);
    }

    return;
}
