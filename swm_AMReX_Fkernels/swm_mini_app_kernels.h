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
    cu(i,j,k) = 0.5*(p(i,j,k) + p(i+1,j,k))*u(i,j,k);
    cv(i,j,k) = 0.5*(p(i,j,k) + p(i,j+1,k))*v(i,j,k);
    z(i,j,k) = (fsdx*(v(i+1,j,k)-v(i,j,k)) + fsdy*(u(i,j+1,k)-u(i,j,k)))/(p(i,j,k)+p(i+1,j,k)+p(i,j+1,k)+p(i+1,j+1,k));
    h(i,j,k) = p(i,j,k) + 0.25*(u(i-1,j,k)*u(i-1,j,k) + u(i,j,k)*u(i,j,k) + v(i,j-1,k)*v(i,j-1,k) + v(i,j,k)*v(i,j,k));
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
    u_new(i,j,k) = u_old(i,j,k) + tdts8 * (z(i,j-1,k)+z(i,j,k)) * (cv(i,j-1,k) + cv(i,j,k) + cv(i+1,j-1,k) + cv(i+1,j,k)) - tdtsdx * (h(i+1,j,k) - h(i,j,k));
    v_new(i,j,k) = v_old(i,j,k) - tdts8 * (z(i-1,j,k)+z(i,j,k)) * (cu(i-1,j,k) + cu(i-1,j+1,k) + cu(i,j,k) + cu(i,j+1,k)) - tdtsdy * (h(i,j+1,k) - h(i,j,k));
    p_new(i,j,k) = p_old(i,j,k) - tdtsdx * (cu(i,j,k) - cu(i-1,j,k)) - tdtsdy * (cv(i,j,k) - cv(i,j-1,k));
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
