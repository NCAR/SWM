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

#endif // SWM_MINI_APP_KERNELS_H_
