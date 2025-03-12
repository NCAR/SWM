module params
  implicit none

  ! Define parameters that correspond to the macros in params.h
  integer, parameter :: M = 512
  integer, parameter :: N = 512
  integer, parameter :: M_LEN = M + 1
  integer, parameter :: N_LEN = N + 1
  integer, parameter :: ITMAX = 4000
  logical, parameter :: L_OUT = .true.
  logical, parameter :: COPY = .false.
  logical, parameter :: VAL_OUT = .false.


end module params

