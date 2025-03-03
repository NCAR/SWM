module params
  implicit none

  ! Define parameters that correspond to the macros in params.h
  integer, parameter :: M = 4096
  integer, parameter :: N = 4096
  integer, parameter :: M_LEN = M + 1
  integer, parameter :: N_LEN = N + 1
  integer, parameter :: ITMAX = 4000
  logical, parameter :: L_OUT = .true.
  logical, parameter :: _COPY_ = .false.
  logical, parameter :: VAL_OUT = .false.


end module params

