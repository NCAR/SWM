module SWM_Fortran_Kernels

  implicit none

contains

  subroutine UpdateIntermediateVariablesKernel(i,j,k,fsdx,fsdy,p,u,v,cu,cv,h,z)

    integer, intent(in) :: i,j,k
    real,    intent(in) :: fsdx, fsdy
    real, dimension(:,:,:), pointer :: p, u, v
    real, dimension(:,:,:), intent(inout) :: cu, cv, h, z

    cu(i+1,j,k) = 0.5 * (p(i+1,j,k) + p(i,j,k)) * u(i+1,j,k)
    cv(i,j+1,k) = 0.5 * (p(i,j+1,k) + p(i,j,k)) * v(i,j+1,k)
    z(i+1,j+1,k) = (fsdx * (v(i+1,j+1,k) - v(i,j+1,k)) - fsdy * (u(i+1,j+1,k) - u(i+1,j,k))) / &
                 (p(i,j,k) + p(i+1,j,k) + p(i+1,j+1,k) + p(i, j+1,k))
    h(i,j,k) = p(i,j,k) + 0.25 * (u(i+1,j,k) * u(i+1,j,k) + u(i,j,k) * u(i,j,k) + &
                                  v(i,j+1,k) * v(i,j+1,k) + v(i,j,k) * v(i,j,k))
  end subroutine UpdateIntermediateVariablesKernel

  subroutine UpdateNewVariablesKernel(i,j,k,tdtsdx,tdtsdy,tdts8,pold,uold,vold,cu,cv,h,z,pnew,unew,vnew)

    integer, intent(in) :: i,j,k
    real, intent(in) :: tdtsdx,tdtsdy,tdts8
    real, dimension(:,:,:), intent(in) :: cu, cv, z, h
    real, dimension(:,:,:), pointer :: pold, uold, vold, pnew, unew, vnew

    unew(i+1,j,k) = uold(i+1,j,k) + &
                  tdts8 * (z(i+1,j+1,k) + z(i+1,j,k)) * (cv(i+1,j+1,k) + cv(i,j+1,k) + cv(i,j,k) + cv(i+1,j,k)) - &
                  tdtsdx * (h(i+1,j,k) - h(i,j,k))
    vnew(i,j+1,k) = vold(i,j+1,k) - &
                  tdts8 * (z(i+1,j+1,k) + z(i,j+1,k)) * (cu(i+1,j+1,k) + cu(i,j+1,k) + cu(i,j,k) + cu(i+1,j,k)) - &
                  tdtsdy * (h(i,j+1,k) - h(i,j,k))
    pnew(i,j,k) = pold(i,j,k) - tdtsdx * (cu(i+1,j,k) - cu(i,j,k)) - tdtsdy * (cv(i,j+1,k) - cv(i,j,k))
  end subroutine UpdateNewVariablesKernel

  subroutine UpdateOldVariablesKernel(i,j,k,alpha,pnew,unew,vnew,p,u,v,pold,uold,vold)

    integer, intent(in) :: i,j,k
    real,    intent(in) :: alpha
    real, dimension(:,:,:), pointer :: pold, uold, vold, p, u, v, pnew, unew, vnew

    real, dimension(size(pold,1),size(pold,2),size(pold,3)) :: pold_cpy, uold_cpy, vold_cpy

    pold_cpy(i,j,k) = pold(i,j,k)
    uold_cpy(i,j,k) = uold(i,j,k)
    vold_cpy(i,j,k) = vold(i,j,k)

    uold(i,j,k) = u(i,j,k) + alpha*(unew(i,j,k) - 2. * u(i,j,k) + uold_cpy(i,j,k))
    vold(i,j,k) = v(i,j,k) + alpha*(vnew(i,j,k) - 2. * v(i,j,k) + vold_cpy(i,j,k))
    pold(i,j,k) = p(i,j,k) + alpha*(pnew(i,j,k) - 2. * p(i,j,k) + pold_cpy(i,j,k))

  end subroutine UpdateOldVariablesKernel

end module SWM_Fortran_Kernels