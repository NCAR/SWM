module shallow_fortran_acc
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public :: fortran_time_loop

contains

    subroutine fortran_time_loop(cu, cv, z, h, u, unew, uold, &
                                v, vnew, vold, p, pnew, pold, &
                                M_LEN, N_LEN, M, N, fsdx, fsdy, &
                                dx, dy, alpha, tdt, time, itmax) bind (C, name="fortran_time_loop")
        
        integer(c_int), intent(in), value :: M_LEN, N_LEN, M, N, itmax
        real(c_double), intent(in), value :: fsdx, fsdy, dx, dy, alpha
        real(c_double), intent(inout) :: tdt, time
        real(c_double), dimension(1:M_LEN,1:N_LEN), intent(inout) :: cu, cv, z, h, &
                                                                    u, unew, uold, &
                                                                    v, vnew, vold, &
                                                                    p, pnew, pold

        ! Local variables
        integer :: i 
        real(c_double) :: tdts8, tdtsdx, tdtsdy
        type(c_ptr) :: u_ptr = c_null_ptr, v_ptr = c_null_ptr, &
                    p_ptr = c_null_ptr, unew_ptr = c_null_ptr, &
                    vnew_ptr = c_null_ptr, pnew_ptr = c_null_ptr

        ! Set up C pointers
        u_ptr = c_loc(u)
        v_ptr = c_loc(v)
        p_ptr = c_loc(p)
        unew_ptr = c_loc(unew)
        vnew_ptr = c_loc(vnew)
        pnew_ptr = c_loc(pnew)

        ! ** Start of time loop ** 
        do i = 1, itmax
            call compute_cu_cv_z_h(M_LEN, N_LEN, M, N, p, u, v, cu, cv, z, h, fsdx, fsdy)
            call apply_bc_cu_cv_z_h(M_LEN, N_LEN, M, N, cu, cv, z, h)
            tdts8 = tdt / 8.d0
            tdtsdx = tdt / dx
            tdtsdy = tdt / dy
            call update_unew_vnew_pnew(M_LEN, N_LEN, M, N, unew, uold, vnew, vold, &
                                    pnew, pold, z, cu, cv, h, tdts8, tdtsdx, tdtsdy)
            call apply_bc_unew_vnew_pnew(M_LEN, N_LEN, M, N, unew, vnew, pnew)
            time = time + tdt
            if ( i > 1 ) then
                call time_smoothing(M_LEN, N_LEN, alpha, u, unew, uold, v, vnew, vold, p, pnew, pold)
            else
                tdt = tdt + tdt
                call first_cycle_copy(M_LEN, N_LEN, u, v, p, uold, vold, pold)
            end if
            !!!!! acc wait for async call !!!!!
            call dswap(u_ptr, unew_ptr)
            call dswap(v_ptr, vnew_ptr)
            call dswap(p_ptr, pnew_ptr)
        end do
        ! ** End of time loop **

        ! Set C pointers to null
        u_ptr = c_null_ptr
        v_ptr = c_null_ptr
        p_ptr = c_null_ptr
        unew_ptr = c_null_ptr
        vnew_ptr = c_null_ptr
        pnew_ptr = c_null_ptr
    end subroutine

    subroutine compute_cu_cv_z_h(M_LEN, N_LEN, M, N, p, u, v, cu, cv, z, h, fsdx, fsdy)
        
        integer(c_int), intent(in) :: M_LEN, N_LEN, M, N
        real(c_double), intent(in) :: fsdx, fsdy
        real(c_double), dimension(1:M_LEN,1:N_LEN), intent(in) :: p, u, v
        real(c_double), dimension(1:M_LEN,1:N_LEN), intent(out) :: cu, cv, z, h

        ! Local variables
        integer :: i, j

        !$acc parallel loop gang vector collapse(2) deviceptr(p,u,v,cu,cv,z,h)
        do j = 1, N
            do i = 1, M
                cu(i+1,j) = 0.5d0 * (p(i+1,j) + p(i,j)) * u(i+1,j)
                cv(i,j+1) = 0.5d0 * (p(i,j+1) + p(i,j)) * v(i,j+1)
                z(i+1,j+1) = (fsdx * (v(i+1,j+1) - v(i,j+1)) - fsdy * (u(i+1,j+1) - u(i+1,j))) / &
                            (p(i,j) + p(i+1,j) + p(i+1,j+1) + p(i,j+1))
                h(i,j) = p(i,j) + 0.25d0 * (u(i+1,j) * u(i+1,j) + u(i,j) * u(i,j) + &
                                            v(i,j+1) * v(i,j+1) + v(i,j) * v(i,j))
            end do
        end do
        !$acc end parallel
    end subroutine

    subroutine apply_bc_cu_cv_z_h(M_LEN, N_LEN, M, N, cu, cv, z, h)

        integer(c_int), intent(in) :: M_LEN, N_LEN, M, N
        real(c_double), dimension(1:M_LEN,1:N_LEN), intent(inout) :: cu, cv, z, h

        ! Local variables
        integer :: i, j

        !$acc parallel loop gang vector deviceptr(cu,cv,z,h)
        do j = 1, N
            cu(1,j) = cu(M_LEN,j)
            cv(M_LEN,j+1) = cv(1,j+1)
            z(1,j+1) = z(M_LEN,j+1)
            h(M_LEN,j) = h(1,j)
        end do
        !$acc end parallel

        !$acc parallel loop gang vector deviceptr(cu,cv,z,h)
        do i = 1, M
            cu(i+1,N_LEN) = cu(i+1,1)
            cv(i,1) = cv(i,N_LEN)
            z(i+1,1) = z(i+1,N_LEN)
            h(i,N_LEN) = h(i,1)
        end do
        !$acc end parallel

        !$acc kernels deviceptr(cu,cv,z,h)
        cu(1,N_LEN) = cu(M_LEN,1)
        cv(M_LEN,1) = cv(1,N_LEN)
        z(1,1) = z(M_LEN,N_LEN)
        h(M_LEN,N_LEN) = h(1,1)
        !$acc end kernels
    end subroutine

    subroutine update_unew_vnew_pnew(M_LEN, N_LEN, M, N, unew, uold, vnew, vold, &
                                    pnew, pold, z, cu, cv, h, tdts8, tdtsdx, tdtsdy)

        integer(c_int), intent(in) :: M_LEN, N_LEN, M, N
        real(c_double), dimension(1:M_LEN,1:N_LEN), intent(out) :: unew, vnew, pnew
        real(c_double), dimension(1:M_LEN,1:N_LEN), intent(in) :: uold, vold, pold, z, cu, cv, h
        real(c_double), intent(in) :: tdts8, tdtsdx, tdtsdy

        ! Local variables
        integer :: i, j

        !$acc parallel loop gang vector collapse(2) deviceptr(unew,vnew,pnew,uold,vold,pold,z,cu,cv,h)
        do j = 1, N
            do i = 1, M
                unew(i+1,j) = uold(i+1,j) + &
                            tdts8 * (z(i+1,j+1) + z(i+1,j)) * (cv(i+1,j+1) + cv(i,j+1) + cv(i,j) + cv(i+1,j)) - &
                            tdtsdx * (h(i+1,j) - h(i,j))
                vnew(i,j+1) = vold(i,j+1) - &
                            tdts8 * (z(i+1,j+1) + z(i,j+1)) * (cu(i+1,j+1) + cu(i,j+1) + cu(i,j) + cu(i+1,j)) - &
                            tdtsdy * (h(i,j+1) - h(i,j))
                pnew(i,j) = pold(i,j) - tdtsdx * (cu(i+1,j) - cu(i,j)) - tdtsdy * (cv(i,j+1) - cv(i,j))
            end do
        end do
        !$acc end parallel
    end subroutine

    subroutine apply_bc_unew_vnew_pnew(M_LEN, N_LEN, M, N, unew, vnew, pnew)

        integer(c_int), intent(in) :: M_LEN, N_LEN, M, N
        real(c_double), dimension(1:M_LEN,1:N_LEN), intent(inout) :: unew, vnew, pnew

        ! Local variables
        integer :: i, j

        !$acc parallel loop gang vector deviceptr(unew, vnew, pnew)
        do j = 1, N
            unew(1,j) = unew(M_LEN,j)
            vnew(M_LEN,j+1) = vnew(1,j+1)
            pnew(M_LEN,j) = pnew(1,j)
        end do
        !$acc end parallel

        !$acc parallel loop gang vector deviceptr(unew, vnew, pnew)
        do i = 1, M
            unew(i+1,N_LEN) = unew(i+1,1)
            vnew(i,1) = vnew(i,N_LEN)
            pnew(i,N_LEN) = pnew(i,1)
        end do
        !$acc end parallel

        !$acc kernels deviceptr(unew, vnew, pnew)
        unew(1,N_LEN) = unew(M_LEN,1)
        vnew(M_LEN,1) = vnew(1,N_LEN)
        pnew(M_LEN,N_LEN) = pnew(1,1)
        !$acc end kernels
    end subroutine

    subroutine time_smoothing(M_LEN, N_LEN, alpha, u, unew, uold, v, vnew, vold, p, pnew, pold)

        integer(c_int), intent(in) :: M_LEN, N_LEN
        real(c_double), intent(in) :: alpha
        real(c_double), dimension(1:M_LEN,1:N_LEN), intent(in) :: u, unew, v, vnew, p, pnew
        real(c_double), dimension(1:M_LEN,1:N_LEN), intent(inout) :: uold, vold, pold

        ! Local variables
        integer :: i, j

        !$acc parallel loop gang vector collapse(2) deviceptr(u,unew,uold,v,vnew,vold,p,pnew,pold)
        do j = 1, N_LEN
            do i = 1, M_LEN
                uold(i,j) = u(i,j) + alpha * (unew(i,j) - 2.d0 * u(i,j) + uold(i,j))
                vold(i,j) = v(i,j) + alpha * (vnew(i,j) - 2.d0 * v(i,j) + vold(i,j))
                pold(i,j) = p(i,j) + alpha * (pnew(i,j) - 2.d0 * p(i,j) + pold(i,j))
            end do
        end do
        !$acc end parallel
    end subroutine

    subroutine first_cycle_copy(M_LEN, N_LEN, u, v, p, uold, vold, pold)

        integer(c_int), intent(in) :: M_LEN, N_LEN
        real(c_double), dimension(1:M_LEN,1:N_LEN), intent(in) :: u, v, p
        real(c_double), dimension(1:M_LEN,1:N_LEN), intent(out) :: uold, vold, pold

        ! Local variables
        integer :: i, j

        !$acc parallel loop gang vector collapse(2) deviceptr(u, v, p, uold, vold, pold)
        do j = 1, N_LEN
            do i = 1, M_LEN
                uold(i,j) = u(i,j)
                vold(i,j) = v(i,j)
                pold(i,j) = p(i,j)
            end do
        end do
        !$acc end parallel
    end subroutine

    subroutine dswap(c_ptr1, c_ptr2)

        type(c_ptr), intent(inout) :: c_ptr1, c_ptr2

        ! Local variable to hold a pointer
        type(c_ptr) :: c_ptr_tmp

        ! Perform the swap of the addresses
        c_ptr_tmp = c_ptr1
        c_ptr1 = c_ptr2
        c_ptr2 = c_ptr_tmp
    end subroutine dswap

end module shallow_fortran_acc