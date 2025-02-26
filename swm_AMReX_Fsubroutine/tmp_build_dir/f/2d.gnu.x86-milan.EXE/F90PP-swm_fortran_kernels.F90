










module SWM_Fortran_Kernels

  implicit none

contains

  subroutine UpdateIntermediateVariablesKernel(p,u,v,cu,cv,h,z,fsdx,fsdy)

    real, dimension(:,:), pointer :: p, u, v
    real, dimension(:,:), intent(inout) :: cu, cv, h, z
    real,    intent(in) :: fsdx, fsdy

    integer :: i,j

    do j=1,size(cu,2)-1
      do i=1,size(cu,1)-1
        cu(i+1,j) = 0.5 * (p(i+1,j) + p(i,j)) * u(i+1,j)
        cv(i,j+1) = 0.5 * (p(i,j+1) + p(i,j)) * v(i,j+1)
        z(i+1,j+1) = (fsdx * (v(i+1,j+1) - v(i,j+1)) - fsdy * (u(i+1,j+1) - u(i+1,j))) / &
                     (p(i,j) + p(i+1,j) + p(i+1,j+1) + p(i, j+1))
        h(i,j) = p(i,j) + 0.25 * (u(i+1,j) * u(i+1,j) + u(i,j) * u(i,j) + &
                                  v(i,j+1) * v(i,j+1) + v(i,j) * v(i,j))
      end do
    end do
  end subroutine UpdateIntermediateVariablesKernel

  subroutine UpdateNewVariablesSub(lo,hi,pold,pold_lo,pold_hi, &
                                            uold,uold_lo,uold_hi, &
                                            vold,vold_lo,vold_hi, &
                                            cu,cu_lo,cu_hi, &
                                            cv,cv_lo,cv_hi, &
                                            h,h_lo,h_hi, &
                                            z,z_lo,z_hi, &
                                            pnew,pnew_lo,pnew_hi, &
                                            unew,unew_lo,unew_hi, &
                                            vnew,vnew_lo,vnew_hi, &
                                            tdtsdx,tdtsdy,tdts8) bind(C,name="UpdateNewVariablesSub")
    use amrex_fort_module, only : amrex_real
    implicit none

    integer, dimension(2) ::  lo, hi
    integer, dimension(2) ::  pold_lo, pold_hi, uold_lo, uold_hi, vold_lo, vold_hi
    integer, dimension(2) ::  cu_lo, cu_hi, cv_lo,cv_hi, h_lo, h_hi, z_lo, z_hi
    integer, dimension(2) ::  pnew_lo, pnew_hi, unew_lo, unew_hi,vnew_lo,vnew_hi
    real(amrex_real), intent(in) :: cu(cu_lo(1):cu_hi(1),cu_lo(2):cu_hi(2))
    real(amrex_real), intent(in) :: cv(cv_lo(1):cv_hi(1),cv_lo(2):cv_hi(2))
    real(amrex_real), intent(in) :: z(z_lo(1):z_hi(1),z_lo(2):z_hi(2))
    real(amrex_real), intent(in) :: h(h_lo(1):h_hi(1),h_lo(2):h_hi(2))
    real :: pold(pold_lo(1):pold_hi(1),pold_lo(2):pold_hi(2))
    real(amrex_real) :: uold(uold_lo(1):uold_hi(1),uold_lo(2):uold_hi(2))
    real(amrex_real) :: vold(vold_lo(1):vold_hi(1),vold_lo(2):vold_hi(2))
    real(amrex_real) :: pnew(pnew_lo(1):pnew_hi(1),pnew_lo(2):pnew_hi(2))
    real(amrex_real) :: unew(unew_lo(1):unew_hi(1),unew_lo(2):unew_hi(2))
    real(amrex_real) :: vnew(vnew_lo(1):vnew_hi(1),vnew_lo(2):vnew_hi(2))

    real(amrex_real), intent(in) :: tdtsdx,tdtsdy,tdts8

    integer :: i,j

!    print *,'UpdateNewVariablesSub: start of subroutine'
    print *,'tdtsdx,tdtsdy,tdts8: ',tdtsdx, tdtsdy, tdts8
!    do j=1,size(unew,2)-1
!      do i=1,size(unew,1)-1
    print *,'pnew(dim=1): ',pnew_lo(1), pnew_hi(2) 
    print *,'pnew(dim=2): ',pnew_lo(2), pnew_hi(2) 
    print *,'j: lo:hi ',lo(2),hi(2)
    print *,'i: lo:hi ',lo(1),hi(1)
    do j=lo(2),hi(2)-1
      do i=lo(1),hi(1)-1
        unew(i+1,j) = uold(i+1,j) + &
                      tdts8 * (z(i+1,j+1) + z(i+1,j)) * (cv(i+1,j+1) + cv(i,j+1) + cv(i,j) + cv(i+1,j)) - &
                      tdtsdx * (h(i+1,j) - h(i,j))
        vnew(i,j+1) = vold(i,j+1) - &
                      tdts8 * (z(i+1,j+1) + z(i,j+1)) * (cu(i+1,j+1) + cu(i,j+1) + cu(i,j) + cu(i+1,j)) - &
                      tdtsdy * (h(i,j+1) - h(i,j))
        pnew(i,j) = pold(i,j) - tdtsdx * (cu(i+1,j) - cu(i,j)) - tdtsdy * (cv(i,j+1) - cv(i,j))
      end do
    end do
!    print *,'UpdateNewVariablesSub: end of subroutine'

  end subroutine UpdateNewVariablesSub



  subroutine UpdateNewVariablesKernel(pold,uold,vold,cu,cv,h,z,pnew,unew,vnew,tdtsdx,tdtsdy,tdts8)

    real, intent(in) :: tdtsdx,tdtsdy,tdts8
    real, dimension(:,:), intent(in) :: cu, cv, z, h
    real, dimension(:,:), pointer :: pold, uold, vold, pnew, unew, vnew

    integer :: i,j

    do j=1,size(unew,2)-1
      do i=1,size(unew,1)-1
        unew(i+1,j) = uold(i+1,j) + &
                      tdts8 * (z(i+1,j+1) + z(i+1,j)) * (cv(i+1,j+1) + cv(i,j+1) + cv(i,j) + cv(i+1,j)) - &
                      tdtsdx * (h(i+1,j) - h(i,j))
        vnew(i,j+1) = vold(i,j+1) - &
                      tdts8 * (z(i+1,j+1) + z(i,j+1)) * (cu(i+1,j+1) + cu(i,j+1) + cu(i,j) + cu(i+1,j)) - &
                      tdtsdy * (h(i,j+1) - h(i,j))
        pnew(i,j) = pold(i,j) - tdtsdx * (cu(i+1,j) - cu(i,j)) - tdtsdy * (cv(i,j+1) - cv(i,j))
      end do
    end do
  end subroutine UpdateNewVariablesKernel

  subroutine UpdateOldVariablesKernel(pnew,unew,vnew,p,u,v,pold,uold,vold,alpha)

    real,    intent(in) :: alpha
    real, dimension(:,:), pointer :: pold, uold, vold, p, u, v, pnew, unew, vnew

    integer :: i,j

    do j=1,size(uold,2)-1
      do i=1,size(uold,1)-1
        uold(i,j) = u(i,j) + alpha*(unew(i,j) - 2. * u(i,j) + uold(i,j))
        vold(i,j) = v(i,j) + alpha*(vnew(i,j) - 2. * v(i,j) + vold(i,j))
        pold(i,j) = p(i,j) + alpha*(pnew(i,j) - 2. * p(i,j) + pold(i,j))
      end do
    end do

  end subroutine UpdateOldVariablesKernel

end module SWM_Fortran_Kernels
