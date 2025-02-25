#define M 256
#define N 256
#define M_LEN (M + 1)
#define N_LEN (N + 1)
#define L_OUT .true.
#define VAL_OUT .false.
#define ITMAX 4000

Program SWM_Fortran_Driver

  use swm_fortran_amrex_kernels, only : UpdateIntermediateVariablesAmrexKernel
  use swm_fortran_amrex_kernels, only : UpdateNewVariablesAmrexKernel
  use swm_fortran_amrex_kernels, only : UpdateOldVariablesAmrexKernel

  implicit none

  ! Solution arrays
  real, dimension(M_LEN,N_LEN,1), target :: u1, u2, u3, &
                                          v1, v2, v3, &
                                          p1, p2, p3

  real, dimension(M_LEN,N_LEN,1) :: cu, &
                                  cv, &
                                  z, &
                                  h, &
                                  psi

  real, dimension(:,:,:), pointer :: u    => NULL(), &
                                     v    => NULL(), &
                                     p    => NULL(), &
                                     unew => NULL(), &
                                     vnew => NULL(), &
                                     pnew => NULL(), &
                                     uold => NULL(), &
                                     vold => NULL(), &
                                     pold => NULL()

  real :: dt, tdt, dx, dy, a, alpha, el, pi
  real :: tpi, di, dj, pcf
  real :: tdts8, tdtsdx, tdtsdy, fsdx, fsdy

  integer, parameter :: mnmin = min(M,N)
  integer :: ncycle
  integer :: i,j

  ! Timer variables
  real :: mfs100, mfs200, mfs300
  real :: t100, t200, t300
  real :: tstart, ctime, tcyc, time, ptime
  real :: c1, c2

  ! Set up pointers
  u    => u1
  v    => v1
  p    => p1
  unew => u2
  vnew => v2
  pnew => p2
  uold => u3
  vold => v3
  pold => p3

  ! Initialization
  dt = 90.
  tdt = dt

  dx = 1.e5
  dy = 1.e5
  fsdx = 4. / dx
  fsdy = 4. / dy

  a = 1.e6
  alpha = 0.001

  el = N * dx
  pi = atan2(0., -1.)
  tpi = pi + pi
  di = tpi / M
  dj = tpi / N
  pcf = pi * pi * a * a / (el * el)

  ! Initial values of the stream function and p
  do j=0,N_LEN-1
    do i=0,M_LEN-1
      psi(i+1,j+1,1) = a * sin((i + .5) * di) * sin((j + .5) * dj)
      p(i+1,j+1,1) = pcf * (cos(2. * (i) * di) + cos(2. * (j) * dj)) + 50000.
    end do
  end do

  ! Initialize velocities
  do j=1,N
    do i=1,M
      u(i+1,j,1) = -(psi(i+1,j+1,1) - psi(i+1,j,1)) / dy
      v(i,j+1,1) = (psi(i+1,j+1,1) - psi(i,j+1,1)) / dx
    end do
  end do

  ! Periodic continuation
  do j=1,N
    u(1,j,1) = u(M_LEN,j,1)
    v(M_LEN,j,1) = v(1,j,1)
  end do

  do i=1,M
    u(i,N_LEN,1) = u(i,1,1)
    v(i,1,1) = v(i,N_LEN,1)
  end do

  u(1,N_LEN,1) = u(M_LEN,1,1)
  v(M_LEN,1,1) = v(1,N_LEN,1)

  do j=1,N_LEN
    do i=1,M_LEN
      uold(i,j,1) = u(i,j,1)
      vold(i,j,1) = v(i,j,1)
      pold(i,j,1) = p(i,j,1)
    end do
  end do

  if ( L_OUT ) then
    write(*,"(A,I0)") " number of points in the x direction ", N
    write(*,"(A,I0)") " number of points in the y direction ", M
    write(*,"(A,F0.6)") " grid spacing in the x direction     ", dx
    write(*,"(A,F0.6)") " grid spacing in the y direction     ", dy
    write(*,"(A,F0.6)") " time step                           ", dt
    write(*,"(A,F0.6)") " time filter parameter               ", alpha

    write(*, "(A)") " initial diagonal elements of p"
    do i=1,mnmin
      write(*, "(F0.6, 1X)", advance="no") p(i,i,1)
    end do

    write(*, "(/,A)") " initial diagonal elements of u"
    do i=1,mnmin
      write(*, "(F0.6, 1X)", advance="no") u(i,i,1)
    end do

    write(*, "(/,A)") " initial diagonal elements of v"
    do i=1,mnmin
      write(*, "(F0.6, 1X)", advance="no") v(i,i,1)
    end do
    write(*,*)
  end if

  ! Start timer
  call cpu_time(tstart)
  time = 0.
  t100 = 0.
  t200 = 0.
  t300 = 0.

  ! Start of time loop
  do ncycle=1,ITMAX

    call cpu_time(c1)
    do j=1,N
      do i=1,M
        call UpdateIntermediateVariablesAmrexKernel(i,j,1,fsdx,fsdy,p,u,v,cu,cv,h,z)
      end do
    end do
    call cpu_time(c2)
    t100 = t100 + (c2 - c1)

    ! Periodic continuation
    do j=1,N
      cu(1,j,1) = cu(M_LEN,j,1)
      cv(M_LEN,j+1,1) = cv(1,j+1,1)
      z(1,j+1,1) = z(M_LEN,j+1,1)
      h(M_LEN,j,1) = h(1,j,1)
    end do

    do i=1,M
      cu(i+1, N_LEN,1) = cu(i+1,1,1)
      cv(i,1,1) = cv(i,N_LEN,1)
      z(i+1,1,1) = z(i+1,N_LEN,1)
      h(i,N_LEN,1) = h(i,1,1)
    end do

    cu(1,N_LEN,1) = cu(M_LEN,1,1)
    cv(M_LEN,1,1) = cv(1,N_LEN,1)
    z(1,1,1) = z(M_LEN,N_LEN,1)
    h(M_LEN,N_LEN,1) = h(1,1,1)

    ! Compute new values of u, v, and p
    tdts8 = tdt / 8.
    tdtsdx = tdt / dx
    tdtsdy = tdt / dy

    call cpu_time(c1)
    do j=1,N
      do i=1,M
        call UpdateNewVariablesAmrexKernel(i,j,1,tdtsdx,tdtsdy,tdts8,pold,uold,vold,cu,cv,h,z,pnew,unew,vnew)
      end do
    end do
    call cpu_time(c2)
    t200 = t200 + (c2-c1)

    ! Periodic continuation
    do j=1,N
      unew(1,j,1) = unew(M_LEN,j,1)
      vnew(M_LEN,j+1,1) = vnew(1,j+1,1)
      pnew(M_LEN,j,1) = pnew(1,j,1)
    end do

    do i=1,M
      unew(i+1,N_LEN,1) = unew(i+1,1,1)
      vnew(i,1,1) = vnew(i,N_LEN,1)
      pnew(i,N_LEN,1) = pnew(i,1,1)
    end do

    unew(1,N_LEN,1) = unew(M_LEN,1,1)
    vnew(M_LEN,1,1) = vnew(1,N_LEN,1)
    pnew(M_LEN,N_LEN,1) = pnew(1,1,1)

    time = time + dt
    if (ncycle > 1) then
      call cpu_time(c1)
      do j=1,N_LEN
        do i=1,M_LEN
          call UpdateOldVariablesAmrexKernel(i,j,1,alpha,pnew,unew,vnew,p,u,v,pold,uold,vold)
        end do
      end do

#ifdef _COPY_
      do j=1,N_LEN
        do i=1,M_LEN
          u(i,j,1) = unew(i,j,1)
          v(i,j,1) = vnew(i,j,1)
          p(i,j,1) = pnew(i,j,1)
        end do
      end do
#else
      call dswap(u, unew)
      call dswap(v, vnew)
      call dswap(p, pnew)
#endif
      call cpu_time(c2)
      t300 = t300 + (c2 - c1)
    else ! ncycle = 1
      tdt = tdt + tdt
      do j=1,N_LEN
        do i=1,N_LEN
          uold(i,j,1) = u(i,j,1)
          vold(i,j,1) = v(i,j,1)
          pold(i,j,1) = p(i,j,1)
        end do
      end do
      call dswap(u, unew)
      call dswap(v, vnew)
      call dswap(p, pnew)
    end if
  end do ! End of time loop

  call dswap(u, unew)
  call dswap(v, vnew)
  call dswap(p, pnew)

  if ( VAL_OUT ) then
    call write_to_file(pnew, 'p.bin')
    call write_to_file(unew, 'u.bin')
    call write_to_file(vnew, 'v.bin')
  end if

  if ( L_OUT ) then
    ptime = time / 3600.

    write(*, "(A,I0,A,F0.6)") " cycle number ", ITMAX, &
                              " model time in hours ", ptime
    write(*, "(A)") " diagonal elements of p"
    do i=1,mnmin
      write(*, "(F0.6,1X)", advance="no") pnew(i,i,1)
    end do
    write(*, "(/,A)") " diagonal elements of u"
    do i=1,mnmin
      write(*, "(F0.6,1X)", advance="no") unew(i,i,1)
    end do
    write(*, "(/,A)") " diagonal elements of v"
    do i=1,mnmin
      write(*, "(F0.6,1X)", advance="no") vnew(i,i,1)
    end do

    mfs100 = 0.
    mfs200 = 0.
    mfs300 = 0.
    ! gdr t100 etc. now an accumulation of all l100 time
    if ( t100 .gt. 0 ) mfs100 = real(ITMAX) * 24. * real(M*N) / t100 / 1000000.
    if ( t200 .gt. 0 ) mfs200 = real(ITMAX) * 26. * real(M*N) / t200 / 1000000.
    if ( t300 .gt. 0 ) mfs300 = real(ITMAX) * 15. * real(M*N) / t300 / 1000000.

    call cpu_time(c2)
    ctime = c2 - tstart
    tcyc = ctime / real(ITMAX)

    write(*, "(/,A,I0,A,F0.6,A,F0.6)") " cycle number ", ITMAX, " total computer time ", ctime, " time per cycle ", tcyc
    write(*, "(A,F0.6,1X,F0.6)") " time and megaflops for loop 100 ", t100, mfs100
    write(*, "(A,F0.6,1X,F0.6)") " time and megaflops for loop 200 ", t200, mfs200
    write(*, "(A,F0.6,1X,F0.6)") " time and megaflops for loop 300 ", t300, mfs300
  end if

contains

  ! This is a copy, not a pointer shuffle
  subroutine dswap(swp1, swp2)

    real, dimension(:,:,:), pointer :: swp1, swp2

    real, dimension(:,:,:), pointer :: swp_tmp

    swp_tmp => swp1
    swp1 => swp2
    swp2 => swp_tmp

  end subroutine dswap

  subroutine write_to_file(array, filename)
    real, dimension(M_LEN,N_LEN,1), intent(in) :: array
    character(len=*),               intent(in) :: filename

    integer :: id

    open(newunit=id, access='stream', status='replace', file=filename)
    ! Write this out in C ordering of array
    do i=1,N_LEN
      do j=1,M_LEN
        write(id) array(i,j,1)
      end do
    end do
    close(id)

  end subroutine write_to_file

End Program SWM_Fortran_Driver