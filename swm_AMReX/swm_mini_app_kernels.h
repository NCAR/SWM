#ifndef SWM_MINI_APP_KERNELS_H_
#define SWM_MINI_APP_KERNELS_H_

#include <AMReX.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void UpdateIntermediateVariablesKernel( const int i, const int j, const int k, const double fsdx, const double fsdy,
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
                               const amrex::Array4<amrex::Real const>& p_old_array,
                               const amrex::Array4<amrex::Real const>& u_old_array,
                               const amrex::Array4<amrex::Real const>& v_old_array,
                               const amrex::Array4<amrex::Real const>& cu_array,
                               const amrex::Array4<amrex::Real const>& cv_array,
                               const amrex::Array4<amrex::Real const>& h_array,
                               const amrex::Array4<amrex::Real const>& z_array,
                               const amrex::Array4<amrex::Real>& p_new_array,
                               const amrex::Array4<amrex::Real>& u_new_array,
                               const amrex::Array4<amrex::Real>& v_new_array)
{
            u_new_array(i,j,k) = u_old_array(i,j,k) + tdts8 * (z_array(i,j-1,k)+z_array(i,j,k)) * (cv_array(i,j-1,k) + cv_array(i,j,k) + cv_array(i+1,j-1,k) + cv_array(i+1,j,k)) - tdtsdx * (h_array(i+1,j,k) - h_array(i,j,k));
            v_new_array(i,j,k) = v_old_array(i,j,k) - tdts8 * (z_array(i-1,j,k)+z_array(i,j,k)) * (cu_array(i-1,j,k) + cu_array(i-1,j+1,k) + cu_array(i,j,k) + cu_array(i,j+1,k)) - tdtsdy * (h_array(i,j+1,k) - h_array(i,j,k));
            p_new_array(i,j,k) = p_old_array(i,j,k) - tdtsdx * (cu_array(i,j,k) - cu_array(i-1,j,k)) - tdtsdy * (cv_array(i,j,k) - cv_array(i,j-1,k));
}

#endif // SWM_MINI_APP_KERNELS_H_
