module SWM_Fortran_Kernels

  implicit none

contains

  subroutine UpdateIntermediateVariablesKernel(fsdx,fsdy,p,u,v,cu,cv,h,z)

    real,    intent(in) :: fsdx, fsdy
    real, dimension(:,:), pointer :: p, u, v
    real, dimension(:,:), intent(inout) :: cu, cv, h, z

    integer :: i,j

    !$acc parallel default(present) 
    !$acc loop collapse(2)
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
    !$acc end parallel

  end subroutine UpdateIntermediateVariablesKernel

  subroutine UpdateNewVariablesKernel(tdtsdx,tdtsdy,tdts8,pold,uold,vold,cu,cv,h,z,pnew,unew,vnew)

    real, intent(in) :: tdtsdx,tdtsdy,tdts8
    real, dimension(:,:), intent(in) :: cu, cv, z, h
    real, dimension(:,:), pointer :: pold, uold, vold, pnew, unew, vnew

    integer :: i,j

    !$acc parallel loop collapse(2) default(present)
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

  subroutine UpdateOldVariablesKernel(alpha,pnew,unew,vnew,p,u,v,pold,uold,vold)

    real,    intent(in) :: alpha
    real, dimension(:,:), pointer :: pold, uold, vold, p, u, v, pnew, unew, vnew

    integer :: i,j

    !$acc parallel loop collapse(2) default(present)
    do j=1,size(uold,2)
      do i=1,size(uold,1)
        uold(i,j) = u(i,j) + alpha*(unew(i,j) - 2. * u(i,j) + uold(i,j))
        vold(i,j) = v(i,j) + alpha*(vnew(i,j) - 2. * v(i,j) + vold(i,j))
        pold(i,j) = p(i,j) + alpha*(pnew(i,j) - 2. * p(i,j) + pold(i,j))
      end do
    end do

  end subroutine UpdateOldVariablesKernel

end module SWM_Fortran_Kernels
