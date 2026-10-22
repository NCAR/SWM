// Kokkos port of module shallow_fortran_acc (OpenACC shallow water time loop).
//
// The arithmetic and the C interface are those of the Fortran version; only the
// kernel structure differs. Per time step the OpenACC code launches 9 kernels
// (cu/cv/z/h, 3 boundary kernels, unew/vnew/pnew, 3 boundary kernels,
// smoothing or first-cycle copy). Here a step is 2 kernels:
//
//   1. cu, cv, z, h: the periodic boundary rows/columns are not copied from
//      the interior afterwards; the thread that computes the interior value
//      also stores it at the boundary location it is copied to.
//   2. unew, vnew, pnew with the same boundary trick, plus the time smoothing
//      (or, on the first step, the copy uold = u) of the same points in the
//      same kernel, since both only need that point's own value.
//
// Every point of each output array is still written exactly once per step.
// The fusion relies on M_LEN == M + 1 and N_LEN == N + 1, which is checked.

#include "shallow_kokkos.h"

#include <Kokkos_Core.hpp>

#include <cstdio>
#include <cstdlib>
#include <utility>

#ifndef SWM_EXEC_SPACE
#define SWM_EXEC_SPACE Kokkos::DefaultExecutionSpace
#endif

// Named namespace (not anonymous): nvcc does not allow extended lambdas in
// functions with internal linkage.
namespace swm {

using ExecSpace = SWM_EXEC_SPACE;
using MemSpace = ExecSpace::memory_space;

// Unmanaged views over the caller's memory. Fortran layout: first index fastest.
// Restrict: the 13 arrays never alias.
using Traits = Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Restrict>;
using Field = Kokkos::View<double **, Kokkos::LayoutLeft, MemSpace, Traits>;
using CField = Kokkos::View<const double **, Kokkos::LayoutLeft, MemSpace, Traits>;

// Iterate::Left on both levels: i (stride 1) is the fastest-varying index, i.e.
// consecutive GPU threads / SIMD lanes touch consecutive addresses.
// Tile sizes are untuned starting points.
constexpr bool kHost = Kokkos::SpaceAccessibility<ExecSpace, Kokkos::HostSpace>::accessible;
constexpr int kTileI = kHost ? 128 : 64;
constexpr int kTileJ = kHost ? 8 : 4;

using Policy = Kokkos::MDRangePolicy<ExecSpace,
                                     Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>,
                                     Kokkos::IndexType<int>>;

// 0-based loop over i = 0..M-1, j = 0..N-1, i.e. Fortran i = 1..M, j = 1..N.
// Fortran X(i+a, j+b) is X(i+a, j+b) here as well (both shifts are -1).
inline Policy interior(int M, int N) { return Policy({0, 0}, {M, N}, {kTileI, kTileJ}); }

KOKKOS_INLINE_FUNCTION double smooth(bool first, double alpha, double cur, double nxt, double old) {
  return first ? cur : cur + alpha * (nxt - 2.0 * cur + old);
}

// Store the new value at (a,b) and the time-smoothed (or first-cycle) old value.
KOKKOS_INLINE_FUNCTION void store(const Field &nw, const Field &old, const CField &cur,
                                  int a, int b, double nx, double ov, bool first, double alpha) {
  nw(a, b) = nx;
  old(a, b) = smooth(first, alpha, cur(a, b), nx, ov);
}

// compute_cu_cv_z_h + apply_bc_cu_cv_z_h
void compute_cu_cv_z_h(int M, int N, CField p, CField u, CField v,
                       Field cu, Field cv, Field z, Field h, double fsdx, double fsdy) {
  Kokkos::parallel_for(
      "swm_cu_cv_z_h", interior(M, N), KOKKOS_LAMBDA(const int i, const int j) {
        const double cu_v = 0.5 * (p(i + 1, j) + p(i, j)) * u(i + 1, j);
        const double cv_v = 0.5 * (p(i, j + 1) + p(i, j)) * v(i, j + 1);
        const double z_v = (fsdx * (v(i + 1, j + 1) - v(i, j + 1)) - fsdy * (u(i + 1, j + 1) - u(i + 1, j))) /
                           (p(i, j) + p(i + 1, j) + p(i + 1, j + 1) + p(i, j + 1));
        const double h_v = p(i, j) + 0.25 * (u(i + 1, j) * u(i + 1, j) + u(i, j) * u(i, j) +
                                             v(i, j + 1) * v(i, j + 1) + v(i, j) * v(i, j));

        // cu(i+1,j): periodic images (0,j) = (M,j), (i+1,N) = (i+1,0), (0,N) = (M,0)
        cu(i + 1, j) = cu_v;
        if (i == M - 1) cu(0, j) = cu_v;
        if (j == 0) cu(i + 1, N) = cu_v;
        if (i == M - 1 && j == 0) cu(0, N) = cu_v;

        // cv(i,j+1): images (M,j+1) = (0,j+1), (i,0) = (i,N), (M,0) = (0,N)
        cv(i, j + 1) = cv_v;
        if (i == 0) cv(M, j + 1) = cv_v;
        if (j == N - 1) cv(i, 0) = cv_v;
        if (i == 0 && j == N - 1) cv(M, 0) = cv_v;

        // z(i+1,j+1): images (0,j+1) = (M,j+1), (i+1,0) = (i+1,N), (0,0) = (M,N)
        z(i + 1, j + 1) = z_v;
        if (i == M - 1) z(0, j + 1) = z_v;
        if (j == N - 1) z(i + 1, 0) = z_v;
        if (i == M - 1 && j == N - 1) z(0, 0) = z_v;

        // h(i,j): images (M,j) = (0,j), (i,N) = (i,0), (M,N) = (0,0)
        h(i, j) = h_v;
        if (i == 0) h(M, j) = h_v;
        if (j == 0) h(i, N) = h_v;
        if (i == 0 && j == 0) h(M, N) = h_v;
      });
}

// update_unew_vnew_pnew + apply_bc_unew_vnew_pnew + time_smoothing / first_cycle_copy
void update_and_smooth(int M, int N, double tdts8, double tdtsdx, double tdtsdy, double alpha, bool first,
                       CField u, CField v, CField p, CField z, CField cu, CField cv, CField h,
                       Field unew, Field vnew, Field pnew, Field uold, Field vold, Field pold) {
  Kokkos::parallel_for(
      "swm_update_smooth", interior(M, N), KOKKOS_LAMBDA(const int i, const int j) {
        const double uo = uold(i + 1, j);
        const double vo = vold(i, j + 1);
        const double po = pold(i, j);

        const double un = uo + tdts8 * (z(i + 1, j + 1) + z(i + 1, j)) *
                                   (cv(i + 1, j + 1) + cv(i, j + 1) + cv(i, j) + cv(i + 1, j)) -
                          tdtsdx * (h(i + 1, j) - h(i, j));
        const double vn = vo - tdts8 * (z(i + 1, j + 1) + z(i, j + 1)) *
                                   (cu(i + 1, j + 1) + cu(i, j + 1) + cu(i, j) + cu(i + 1, j)) -
                          tdtsdy * (h(i, j + 1) - h(i, j));
        const double pn = po - tdtsdx * (cu(i + 1, j) - cu(i, j)) - tdtsdy * (cv(i, j + 1) - cv(i, j));

        // u at (i+1,j): images (0,j), (i+1,N), (0,N)
        store(unew, uold, u, i + 1, j, un, uo, first, alpha);
        if (i == M - 1) store(unew, uold, u, 0, j, un, uold(0, j), first, alpha);
        if (j == 0) store(unew, uold, u, i + 1, N, un, uold(i + 1, N), first, alpha);
        if (i == M - 1 && j == 0) store(unew, uold, u, 0, N, un, uold(0, N), first, alpha);

        // v at (i,j+1): images (M,j+1), (i,0), (M,0)
        store(vnew, vold, v, i, j + 1, vn, vo, first, alpha);
        if (i == 0) store(vnew, vold, v, M, j + 1, vn, vold(M, j + 1), first, alpha);
        if (j == N - 1) store(vnew, vold, v, i, 0, vn, vold(i, 0), first, alpha);
        if (i == 0 && j == N - 1) store(vnew, vold, v, M, 0, vn, vold(M, 0), first, alpha);

        // p at (i,j): images (M,j), (i,N), (M,N)
        store(pnew, pold, p, i, j, pn, po, first, alpha);
        if (i == 0) store(pnew, pold, p, M, j, pn, pold(M, j), first, alpha);
        if (j == 0) store(pnew, pold, p, i, N, pn, pold(i, N), first, alpha);
        if (i == 0 && j == 0) store(pnew, pold, p, M, N, pn, pold(M, N), first, alpha);
      });
}

void ensure_kokkos() {
  if (Kokkos::is_finalized()) {
    std::fprintf(stderr, "kokkos_time_loop: Kokkos has already been finalized\n");
    std::abort();
  }
  if (!Kokkos::is_initialized()) {
    Kokkos::initialize();
    std::atexit(kokkos_swm_finalize);
  }
}

}  // namespace swm

extern "C" {

void kokkos_swm_init(void) {
  if (!Kokkos::is_initialized()) Kokkos::initialize();
}

void kokkos_swm_finalize(void) {
  if (Kokkos::is_initialized() && !Kokkos::is_finalized()) Kokkos::finalize();
}

void kokkos_time_loop(double *cu, double *cv, double *z, double *h,
                      double *u, double *unew, double *uold,
                      double *v, double *vnew, double *vold,
                      double *p, double *pnew, double *pold,
                      int M_LEN, int N_LEN, int M, int N,
                      double fsdx, double fsdy, double dx, double dy,
                      double alpha, double *tdt, double *sim_time, int itmax) {
  if (M_LEN != M + 1 || N_LEN != N + 1 || M < 1 || N < 1) {
    std::fprintf(stderr,
                 "kokkos_time_loop: requires M_LEN == M+1, N_LEN == N+1, M,N >= 1 "
                 "(got M_LEN=%d N_LEN=%d M=%d N=%d)\n",
                 M_LEN, N_LEN, M, N);
    std::abort();
  }
  swm::ensure_kokkos();

  using swm::Field;
  Field f_cu(cu, M_LEN, N_LEN), f_cv(cv, M_LEN, N_LEN), f_z(z, M_LEN, N_LEN), f_h(h, M_LEN, N_LEN);
  Field f_u(u, M_LEN, N_LEN), f_unew(unew, M_LEN, N_LEN), f_uold(uold, M_LEN, N_LEN);
  Field f_v(v, M_LEN, N_LEN), f_vnew(vnew, M_LEN, N_LEN), f_vold(vold, M_LEN, N_LEN);
  Field f_p(p, M_LEN, N_LEN), f_pnew(pnew, M_LEN, N_LEN), f_pold(pold, M_LEN, N_LEN);

  for (int it = 1; it <= itmax; ++it) {
    swm::compute_cu_cv_z_h(M, N, f_p, f_u, f_v, f_cu, f_cv, f_z, f_h, fsdx, fsdy);

    const double tdts8 = *tdt / 8.0;
    const double tdtsdx = *tdt / dx;
    const double tdtsdy = *tdt / dy;

    swm::update_and_smooth(M, N, tdts8, tdtsdx, tdtsdy, alpha, it == 1,
                           f_u, f_v, f_p, f_z, f_cu, f_cv, f_h,
                           f_unew, f_vnew, f_pnew, f_uold, f_vold, f_pold);

    *sim_time += *tdt;
    if (it == 1) *tdt += *tdt;  // first cycle: forward step, then double the time step

    // Pointer swap only; no host/device synchronization needed between steps.
    std::swap(f_u, f_unew);
    std::swap(f_v, f_vnew);
    std::swap(f_p, f_pnew);
  }

  Kokkos::fence("kokkos_time_loop: end of time loop");
}

}  // extern "C"
