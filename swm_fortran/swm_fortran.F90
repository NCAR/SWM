#define M 256
#define N 256
#define M_LEN (M + 1)
#define N_LEN (N + 1)
#define L_OUT .true.
#define ITMAX 4000

Program SWM_Fortran

  implicit none

  ! Solution arrays
  real, dimension(:,:), pointer :: u => NULL(), &
                                   v => NULL(), &
                                   p => NULL(), &
                                   unew => NULL(), &
                                   vnew => NULL(), &
                                   pnew => NULL(), &
                                   uold => NULL(), &
                                   vold => NULL(), &
                                   pold => NULL(), &
                                   cu => NULL(), &
                                   cv => NULL(), &
                                   z => NULL(), &
                                   h => NULL(), &
                                   psi => NULL()

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

  ! Allocate memory
  allocate(u(M_LEN,N_LEN), &
           v(M_LEN,N_LEN), &
           p(M_LEN,N_LEN), &
           unew(M_LEN,N_LEN), &
           vnew(M_LEN,N_LEN), &
           pnew(M_LEN,N_LEN), &
           uold(M_LEN,N_LEN), &
           vold(M_LEN,N_LEN), &
           pold(M_LEN,N_LEN), &
           cu(M_LEN,N_LEN), &
           cv(M_LEN,N_LEN), &
           z(M_LEN,N_LEN), &
           h(M_LEN,N_LEN), &
           psi(M_LEN,N_LEN), &
           source=0.)

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
      psi(i+1,j+1) = a * sin((i + .5) * di) * sin((j + .5) * dj)
      p(i+1,j+1) = pcf * (cos(2. * (i) * di) + cos(2. * (j) * dj)) + 50000.
    end do
  end do

  ! Initialize velocities
  do j=1,N
    do i=1,M
      u(i+1,j) = -(psi(i+1,j+1) - psi(i+1,j)) / dy
      v(i,j+1) = (psi(i+1,j+1) - psi(i,j+1)) / dx
    end do
  end do

  ! Periodic continuation
  do j=1,N
    u(1,j) = u(M_LEN,j)
    v(M_LEN,j) = v(1,j)
  end do

  do i=1,M
    u(i,N_LEN) = u(i,1)
    v(i,1) = v(i,N_LEN)
  end do

  u(1,N_LEN) = u(M_LEN,1)
  v(M_LEN,1) = v(1,N_LEN)

  do j=1,N_LEN
    do i=1,M_LEN
      uold(i,j) = u(i,j)
      vold(i,j) = v(i,j)
      pold(i,j) = p(i,j)
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
      write(*, "(F0.6, 1X)", advance="no") p(i,i)
    end do

    write(*, "(/,A)") " initial diagonal elements of u"
    do i=1,mnmin
      write(*, "(F0.6, 1X)", advance="no") u(i,i)
    end do

    write(*, "(/,A)") " initial diagonal elements of v"
    do i=1,mnmin
      write(*, "(F0.6, 1X)", advance="no") v(i,i)
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
        cu(i+1,j) = 0.5 * (p(i+1,j) + p(i,j)) * u(i+1,j)
        cv(i,j+1) = 0.5 * (p(i,j+1) + p(i,j)) * v(i,j+1)
        z(i+1,j+1) = (fsdx * (v(i+1,j+1) - v(i,j+1)) - fsdy * (u(i+1,j+1) - u(i+1,j))) / &
                     (p(i,j) + p(i+1,j) + p(i+1,j+1) + p(i, j+1))
        h(i,j) = p(i,j) + 0.25 * (u(i+1,j) * u(i+1,j) + u(i,j) * u(i,j) + &
                                  v(i,j+1) * v(i,j+1) + v(i,j) * v(i,j))
      end do
    end do
    call cpu_time(c2)
    t100 = t100 + (c2 - c1)

    ! Periodic continuation
    do j=1,N
      cu(1,j) = cu(M_LEN,j)
      cv(M_LEN,j+1) = cv(1,j+1)
      z(1,j+1) = z(M_LEN,j+1)
      h(M_LEN,j) = h(1,j)
    end do

    do i=1,M
      cu(i+1, N_LEN) = cu(i+1,1)
      cv(i,1) = cv(i,N_LEN)
      z(i+1,1) = z(i+1,N_LEN)
      h(i,N_LEN) = h(i,1)
    end do

    cu(1,N_LEN) = cu(M_LEN,1)
    cv(M_LEN,1) = cv(1,N_LEN)
    z(1,1) = z(M_LEN,N_LEN)
    h(M_LEN,N_LEN) = h(1,1)

    ! Compute new values of u, v, and p
    tdts8 = tdt / 8.
    tdtsdx = tdt / dx
    tdtsdy = tdt / dy

    call cpu_time(c1)
    do j=1,N
      do i=1,M
        unew(i+1,j) = uold(i+1,j) + &
                      tdts8 * (z(i+1,j+1) + z(i+1,j)) * (cv(i+1,j+1) + cv(i,j+1) + cv(i,j) + cv(i+1,j)) - &
                      tdtsdx * (h(i+1,j) - h(i,j))
        vnew(i,j+1) = vold(i,j+1) - &
                      tdts8 * (z(i+1,j+1) + z(i,j+1)) * (cu(i+1,j+1) + cu(i,j+1) + cu(i,j) + cu(i+1,j)) - &
                      tdtsdy * (h(i,j+1) - h(i,j))
        pnew(i,j) = pold(i,j) - tdtsdx * (cu(i+1,j) - cu(i,j)) - tdtsdy * (cv(i,j+1) - cv(i,j))
      end do
    end do
    call cpu_time(c2)
    t200 = t200 + (c2-c1)

    ! Periodic continuation
    do j=1,N
      unew(1,j) = unew(M_LEN,j)
      vnew(M_LEN,j+1) = vnew(1,j+1)
      pnew(M_LEN,j) = pnew(1,j)
    end do

    do i=1,M
      unew(i+1,N_LEN) = unew(i+1,1)
      vnew(i,1) = vnew(i,N_LEN)
      pnew(i,N_LEN) = pnew(i,1)
    end do

    unew(1,N_LEN) = unew(M_LEN,1)
    vnew(M_LEN,1) = vnew(1,N_LEN)
    pnew(M_LEN,N_LEN) = pnew(1,1)

    time = time + dt
    if (ncycle > 1) then
      call cpu_time(c1)
      do j=1,N_LEN
        do i=1,M_LEN
          uold(i,j) = u(i,j) + alpha*(unew(i,j) - 2. * u(i,j) + uold(i,j))
          vold(i,j) = v(i,j) + alpha*(vnew(i,j) - 2. * v(i,j) + vold(i,j))
          pold(i,j) = p(i,j) + alpha*(pnew(i,j) - 2. * p(i,j) + pold(i,j))
        end do
      end do

#ifdef _COPY_
      do j=1,N_LEN
        do i=1,M_LEN
          u(i,j) = unew(i,j)
          v(i,j) = vnew(i,j)
          p(i,j) = pnew(i,j)
        end do
      end do
#else
      call dswap(u, unew)
      call dswap(v, vnew)
      call dswap(p, pnew)
#endif
      call cpu_time(c2)
      t300 = t300 + (c2 - c1)
    else
      tdt = tdt + tdt
      do j=1,N_LEN
        do i=1,N_LEN
          uold(i,j) = u(i,j)
          vold(i,j) = v(i,j)
          pold(i,j) = p(i,j)
        end do
      end do
      call dswap(u, unew)
      call dswap(v, vnew)
      call dswap(p, pnew)
    end if
  end do

  call dswap(u, unew)
  call dswap(v, vnew)
  call dswap(p, pnew)

  if ( L_OUT ) then
    ptime = time / 3600.

    write(*, "(A,I0,A,F0.6)") " cycle number ", ITMAX, &
                              " model time in hours ", ptime
    write(*, "(A)") " diagonal elements of p"
    do i=1,mnmin
      write(*, "(F0.6,1X)", advance="no") pnew(i,i)
    end do
    write(*, "(/,A)") " diagonal elements of u"
    do i=1,mnmin
      write(*, "(F0.6,1X)", advance="no") unew(i,i)
    end do
    write(*, "(/,A)") " diagonal elements of v"
    do i=1,mnmin
      write(*, "(F0.6,1X)", advance="no") vnew(i,i)
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

  deallocate(u, v, p, unew, vnew, pnew, uold, vold, pold, cu, cv, z, h, psi)

contains

  subroutine dswap(a, b)

    real, dimension(:,:), pointer :: a, b

    real, dimension(:,:), pointer :: c

    c => a
    a => b
    b => c

  end subroutine dswap

End Program SWM_Fortran