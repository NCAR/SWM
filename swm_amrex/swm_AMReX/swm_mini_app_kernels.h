#ifndef SWM_MINI_APP_KERNELS_H_
#define SWM_MINI_APP_KERNELS_H_

#include <AMReX.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void UpdateIntermediateVariablesKernel( const int i, const int j, const int k,
                                        const double fsdx, const double fsdy,
                                        const amrex::Array4<amrex::Real const>& p,
                                        const amrex::Array4<amrex::Real const>& u,
                                        const amrex::Array4<amrex::Real const>& v,
                                        const amrex::Array4<amrex::Real>& cu,
                                        const amrex::Array4<amrex::Real>& cv,
                                        const amrex::Array4<amrex::Real>& h,
                                        const amrex::Array4<amrex::Real>& z)
{
    amrex::Real p_ijk = p(i,j,k);
    amrex::Real p_i1jk = p(i+1,j,k);
    amrex::Real p_ij1k = p(i,j+1,k);
    amrex::Real p_i1j1k = p(i+1,j+1,k);
    amrex::Real u_ijk = u(i,j,k);
    amrex::Real u_im1jk = u(i-1,j,k);
    amrex::Real u_ij1k = u(i,j+1,k);
    amrex::Real v_ijk = v(i,j,k);
    amrex::Real v_i1jk = v(i+1,j,k);
    amrex::Real v_ijm1k = v(i,j-1,k);
    cu(i,j,k) = 0.5*(p_ijk + p_i1jk)*u_ijk;
    cv(i,j,k) = 0.5*(p_ijk + p_ij1k)*v_ijk;
    z(i,j,k) = (fsdx*(v_i1jk - v_ijk) - fsdy*(u_ij1k - u_ijk))/(p_ijk + p_i1jk + p_ij1k + p_i1j1k);
    h(i,j,k) = p_ijk + 0.25*(u_im1jk*u_im1jk + u_ijk*u_ijk + v_ijm1k*v_ijm1k + v_ijk*v_ijk);
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void UpdateNewVariablesKernel( const int i, const int j, const int k, 
                               const double tdtsdx, const double tdtsdy, const double tdts8,
                               const amrex::Array4<amrex::Real const>& p_old,
                               const amrex::Array4<amrex::Real const>& u_old,
                               const amrex::Array4<amrex::Real const>& v_old,
                               const amrex::Array4<amrex::Real const>& cu,
                               const amrex::Array4<amrex::Real const>& cv,
                               const amrex::Array4<amrex::Real const>& h,
                               const amrex::Array4<amrex::Real const>& z,
                               const amrex::Array4<amrex::Real>& p_new,
                               const amrex::Array4<amrex::Real>& u_new,
                               const amrex::Array4<amrex::Real>& v_new)
{
    amrex::Real z_ijk = z(i,j,k);
    amrex::Real z_ijm1k = z(i,j-1,k);
    amrex::Real z_im1jk = z(i-1,j,k);
    amrex::Real cv_ijk = cv(i,j,k);
    amrex::Real cv_ijm1k = cv(i,j-1,k);
    amrex::Real cv_i1jk = cv(i+1,j,k);
    amrex::Real cv_i1jm1k = cv(i+1,j-1,k);
    amrex::Real cu_ijk = cu(i,j,k);
    amrex::Real cu_im1jk = cu(i-1,j,k);
    amrex::Real cu_ij1k = cu(i,j+1,k);
    amrex::Real cu_im1j1k = cu(i-1,j+1,k);
    amrex::Real h_ijk = h(i,j,k);
    amrex::Real h_i1jk = h(i+1,j,k);
    amrex::Real h_ij1k = h(i,j+1,k);
    u_new(i,j,k) = u_old(i,j,k) + tdts8 * (z_ijm1k + z_ijk) * (cv_ijm1k + cv_ijk + cv_i1jm1k + cv_i1jk) - tdtsdx * (h_i1jk - h_ijk);
    v_new(i,j,k) = v_old(i,j,k) - tdts8 * (z_im1jk + z_ijk) * (cu_im1jk + cu_im1j1k + cu_ijk + cu_ij1k) - tdtsdy * (h_ij1k - h_ijk);
    p_new(i,j,k) = p_old(i,j,k) - tdtsdx * (cu_ijk - cu_im1jk) - tdtsdy * (cv_ijk - cv_ijm1k);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void UpdateOldVariablesKernel( const int i, const int j, const int k, 
                               const double alpha,
                               const amrex::Array4<amrex::Real const>& p,
                               const amrex::Array4<amrex::Real const>& u,
                               const amrex::Array4<amrex::Real const>& v,
                               const amrex::Array4<amrex::Real const>& p_new,
                               const amrex::Array4<amrex::Real const>& u_new,
                               const amrex::Array4<amrex::Real const>& v_new,
                               const amrex::Array4<amrex::Real>& p_old,
                               const amrex::Array4<amrex::Real>& u_old,
                               const amrex::Array4<amrex::Real>& v_old)
{
    amrex::Real u_old_temp = u_old(i,j,k);
    amrex::Real v_old_temp = v_old(i,j,k);
    amrex::Real p_old_temp = p_old(i,j,k);

    u_old(i,j,k) = u(i,j,k) + alpha * (u_new(i,j,k) - 2.0*u(i,j,k) + u_old_temp);
    v_old(i,j,k) = v(i,j,k) + alpha * (v_new(i,j,k) - 2.0*v(i,j,k) + v_old_temp);
    p_old(i,j,k) = p(i,j,k) + alpha * (p_new(i,j,k) - 2.0*p(i,j,k) + p_old_temp);
}

#endif // SWM_MINI_APP_KERNELS_H_
