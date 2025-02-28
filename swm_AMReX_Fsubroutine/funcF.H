#ifndef FUNCF_H_
#define FUNCF_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    // Fortran subroutines

    void UpdateNewVariablesSub(const int* lo, const int* hi,
           const amrex::Real* pold, const int* pold_lo, const int* pold_hi,
           const amrex::Real* uold, const int* uold_lo, const int* uold_hi,
           const amrex::Real* vold, const int* vold_lo, const int* vold_hi,
           const amrex::Real* cu, const int* cu_lo, const int* cu_hi,
           const amrex::Real* cv, const int* cv_lo, const int* cv_hi,
           const amrex::Real* h, const int*   h_lo, const int*  h_hi,
           const amrex::Real* z, const int*   z_lo, const int*  z_hi,
           amrex::Real* pnew, const int* pnew_lo, const int* pnew_hi,
           amrex::Real* unew, const int* unew_lo, const int* unew_hi,
           amrex::Real* vnew, const int* vnew_lo, const int* vnew_hi,
           const amrex::Real& tdtsdx, const amrex::Real& tdtsdy, const amrex::Real& tdts8);

    void UpdateIntermediateVariablesSub(const int* lo, const int* hi,
           const amrex::Real* p,  const int* p_lo, const int* p_hi,
           const amrex::Real* u,  const int* u_lo, const int* u_hi,
           const amrex::Real* v,  const int* v_lo, const int* v_hi,
           amrex::Real* cu, const int* cu_lo, const int* cu_hi,
           amrex::Real* cv, const int* cv_lo, const int* cv_hi,
           amrex::Real* h,  const int*  h_lo, const int*  h_hi,
           amrex::Real* z,  const int*  z_lo, const int*  z_hi,
           const amrex::Real& fsdx, const amrex::Real& fsdy);

    void UpdateOldVariablesSub(const int* lo, const int* hi,
           const amrex::Real* pnew, const int* pnew_lo, const int* pnew_hi,
           const amrex::Real* unew, const int* unew_lo, const int* unew_hi,
           const amrex::Real* vnew, const int* vnew_lo, const int* vnew_hi,
           const amrex::Real* p, const int* p_lo, const int* p_hi,
           const amrex::Real* u, const int* u_lo, const int* u_hi,
           const amrex::Real* v, const int* v_lo, const int* v_hi,
           amrex::Real* pold, const int* pold_lo, const int* pold_hi,
           amrex::Real* uold, const int* uold_lo, const int* uold_hi,
           amrex::Real* vold, const int* vold_lo, const int* vold_hi,
           const amrex::Real& alpha);

#ifdef __cplusplus
}
#endif

#endif
