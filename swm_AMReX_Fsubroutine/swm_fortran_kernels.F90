module SWM_Fortran_Kernels

  implicit none

contains

  subroutine UpdateIntermediateVariablesSub(lo,hi, &
               p,   p_lo,  p_hi,  u,  u_lo,  u_hi, v, v_lo, v_hi, &
               cu, cu_lo, cu_hi, cv, cv_lo, cv_hi, h, h_lo, h_hi, &
               z,   z_lo,  z_hi, fsdx, fsdy) bind(C,name="UpdateIntermediateVariablesSub")

    use amrex_fort_module, only : amrex_real
    implicit none

    integer, dimension(2) ::  lo, hi
    integer, dimension(2) ::  p_lo, p_hi, u_lo, u_hi, v_lo, v_hi
    integer, dimension(2) ::  cu_lo, cu_hi, cv_lo,cv_hi, h_lo, h_hi, z_lo, z_hi
    real(amrex_real), intent(out) :: cu(cu_lo(1):cu_hi(1),cu_lo(2):cu_hi(2))
    real(amrex_real), intent(out) :: cv(cv_lo(1):cv_hi(1),cv_lo(2):cv_hi(2))
    real(amrex_real), intent(out) :: z(z_lo(1):z_hi(1),z_lo(2):z_hi(2))
    real(amrex_real), intent(out) :: h(h_lo(1):h_hi(1),h_lo(2):h_hi(2))
    real(amrex_real), intent(in) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2))
    real(amrex_real), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2))
    real(amrex_real), intent(in) :: v(v_lo(1):v_hi(1),v_lo(2):v_hi(2))
    real(amrex_real), intent(in)    :: fsdx, fsdy

    integer :: i,j

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        cu(i,j) = 0.5 * (p(i,j) + p(i+1,j)) * u(i,j)
        cv(i,j) = 0.5 * (p(i,j) + p(i,j+1)) * v(i,j)
        z(i,j)  = (fsdx * (v(i+1,j) - v(i,j)) - fsdy * (u(i,j+1) - u(i,j))) / &
                     (p(i,j) + p(i+1,j) + p(i,j+1) + p(i+1, j+1))
        h(i,j)  = p(i,j) + 0.25 * (u(i-1,j) * u(i-1,j) + u(i,j) * u(i,j) + &
                                  v(i,j-1) * v(i,j-1) + v(i,j) * v(i,j))
      end do
    end do
  end subroutine UpdateIntermediateVariablesSub

  subroutine UpdateNewVariablesSub(lo,hi, &
               pold, pold_lo, pold_hi, uold, uold_lo, uold_hi, vold, vold_lo, vold_hi, &
                 cu,   cu_lo,   cu_hi,   cv,   cv_lo,   cv_hi,    h,    h_lo,    h_hi, &
                  z,    z_lo,    z_hi, pnew, pnew_lo, pnew_hi, unew, unew_lo,  unew_hi, &
               vnew, vnew_lo,  vnew_hi, tdtsdx,tdtsdy,tdts8) bind(C,name="UpdateNewVariablesSub")

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
    real(amrex_real), intent(in) :: pold(pold_lo(1):pold_hi(1),pold_lo(2):pold_hi(2))
    real(amrex_real), intent(in) :: uold(uold_lo(1):uold_hi(1),uold_lo(2):uold_hi(2))
    real(amrex_real), intent(in) :: vold(vold_lo(1):vold_hi(1),vold_lo(2):vold_hi(2))
    real(amrex_real) :: pnew(pnew_lo(1):pnew_hi(1),pnew_lo(2):pnew_hi(2))
    real(amrex_real) :: unew(unew_lo(1):unew_hi(1),unew_lo(2):unew_hi(2))
    real(amrex_real) :: vnew(vnew_lo(1):vnew_hi(1),vnew_lo(2):vnew_hi(2))
    real(amrex_real), intent(in) :: tdtsdx,tdtsdy,tdts8

    integer :: i,j

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        unew(i,j) = uold(i,j) + &
                      tdts8 * (z(i,j-1) + z(i,j)) * (cv(i,j-1) + cv(i,j) + cv(i+1,j-1) + cv(i+1,j)) - &
                      tdtsdx * (h(i+1,j) - h(i,j))
        vnew(i,j) = vold(i,j) - &
                      tdts8 * (z(i-1,j) + z(i,j)) * (cu(i-1,j) + cu(i-1,j+1) + cu(i,j) + cu(i,j+1)) - &
                      tdtsdy * (h(i,j+1) - h(i,j))
        pnew(i,j) = pold(i,j) - tdtsdx * (cu(i,j) - cu(i-1,j)) - tdtsdy * (cv(i,j) - cv(i,j-1))
      end do
    end do

  end subroutine UpdateNewVariablesSub

  subroutine UpdateOldVariablesSub(lo,hi, &
               pnew, pnew_lo, pnew_hi, unew, unew_lo, unew_hi, vnew, vnew_lo, vnew_hi, &
                  p,    p_lo,    p_hi,    u,    u_lo,    u_hi,    v,    v_lo,    v_hi, &
               pold, pold_lo, pold_hi, uold, uold_lo, uold_hi, vold, vold_lo, vold_hi, &
               alpha) bind(C,name="UpdateOldVariablesSub")

    use amrex_fort_module, only : amrex_real
    implicit none

    integer, dimension(2) ::  lo, hi
    integer, dimension(2) ::  pold_lo, pold_hi, uold_lo, uold_hi, vold_lo, vold_hi
    integer, dimension(2) ::  pnew_lo, pnew_hi, unew_lo, unew_hi, vnew_lo, vnew_hi
    integer, dimension(2) ::  p_lo, p_hi, u_lo, u_hi, v_lo, v_hi

    real(amrex_real), intent(in) :: pnew(pnew_lo(1):pnew_hi(1),pnew_lo(2):pnew_hi(2))
    real(amrex_real), intent(in) :: unew(unew_lo(1):unew_hi(1),unew_lo(2):unew_hi(2))
    real(amrex_real), intent(in) :: vnew(vnew_lo(1):vnew_hi(1),vnew_lo(2):vnew_hi(2))

    real(amrex_real), intent(in)    :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2))
    real(amrex_real), intent(in)    :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2))
    real(amrex_real), intent(in)    :: v(v_lo(1):v_hi(1),v_lo(2):v_hi(2))

    real(amrex_real) :: pold(pold_lo(1):pold_hi(1),pold_lo(2):pold_hi(2))
    real(amrex_real) :: uold(uold_lo(1):uold_hi(1),uold_lo(2):uold_hi(2))
    real(amrex_real) :: vold(vold_lo(1):vold_hi(1),vold_lo(2):vold_hi(2))
    real(amrex_real), intent(in) :: alpha

    integer :: i,j

    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        uold(i,j) = u(i,j) + alpha*(unew(i,j) - 2. * u(i,j) + uold(i,j))
        vold(i,j) = v(i,j) + alpha*(vnew(i,j) - 2. * v(i,j) + vold(i,j))
        pold(i,j) = p(i,j) + alpha*(pnew(i,j) - 2. * p(i,j) + pold(i,j))
      end do
    end do

  end subroutine UpdateOldVariablesSub

end module SWM_Fortran_Kernels
