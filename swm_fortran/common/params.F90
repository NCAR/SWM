module params
  implicit none

  ! Define parameters that correspond to the macros in params.h
  integer, parameter :: M = MNUM
  integer, parameter :: N = NNUM
  integer, parameter :: M_LEN = M + 1
  integer, parameter :: N_LEN = N + 1
  integer, parameter :: ITMAX = 4000
  logical, parameter :: L_OUT = .false.
  logical, parameter :: COPY = .false.
  logical, parameter :: VAL_OUT = .false.


end module params

