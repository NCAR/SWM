/* C interface to the Kokkos port of shallow_fortran_acc::fortran_time_loop.
 *
 * Same argument list, order and by-value / by-reference conventions as the
 * Fortran bind(C) routine, so a C/C++ driver only needs to change the name of
 * the function it calls.
 *
 * Arrays are M_LEN x N_LEN doubles, column-major (Fortran order, first index
 * contiguous). The pointers must be addressable from Kokkos' default memory
 * space: device pointers (cudaMalloc / acc_malloc / Kokkos allocations) when
 * Kokkos is built with CUDA/HIP/SYCL, ordinary host pointers for Serial/OpenMP.
 *
 * Requires M_LEN == M + 1 and N_LEN == N + 1.
 *
 * Like the Fortran routine, the u/unew, v/vnew and p/pnew buffers are swapped
 * on every step but only inside the routine, so on return the latest fields
 * are in u, v, p when itmax is even and in unew, vnew, pnew when itmax is odd.
 * tdt and sim_time are updated in place.
 */
#ifndef SHALLOW_KOKKOS_H
#define SHALLOW_KOKKOS_H

#ifdef __cplusplus
extern "C" {
#endif

void kokkos_time_loop(double *cu, double *cv, double *z, double *h,
                      double *u, double *unew, double *uold,
                      double *v, double *vnew, double *vold,
                      double *p, double *pnew, double *pold,
                      int M_LEN, int N_LEN, int M, int N,
                      double fsdx, double fsdy, double dx, double dy,
                      double alpha, double *tdt, double *sim_time, int itmax);

/* Optional. kokkos_time_loop() initializes Kokkos on first use (and finalizes
 * it at exit) if the driver has not already called Kokkos::initialize(). */
void kokkos_swm_init(void);
void kokkos_swm_finalize(void);

#ifdef __cplusplus
}
#endif

#endif /* SHALLOW_KOKKOS_H */
