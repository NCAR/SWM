// Checks kokkos_time_loop() against a literal, unfused host translation of the
// Fortran/OpenACC module (separate boundary loops, separate smoothing pass,
// 1-based column-major indexing). All 13 arrays are compared, so a missed or
// misplaced boundary point shows up.
//
// Fields are random-ish at every index including the boundaries, so the test
// does not rely on the input being periodic.

#include "shallow_kokkos.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#ifndef SWM_EXEC_SPACE
#define SWM_EXEC_SPACE Kokkos::DefaultExecutionSpace
#endif
using MemSpace = SWM_EXEC_SPACE::memory_space;
using DevView = Kokkos::View<double **, Kokkos::LayoutLeft, MemSpace>;
using HostView = Kokkos::View<double **, Kokkos::LayoutLeft, Kokkos::HostSpace>;

static const char *kNames[13] = {"cu", "cv", "z", "h", "u", "unew", "uold",
                                 "v", "vnew", "vold", "p", "pnew", "pold"};
enum { CU, CV, Z, H, U, UNEW, UOLD, V, VNEW, VOLD, P, PNEW, POLD };

// Fortran-style access: 1-based, first index contiguous.
#define A(f, i, j) f[((i)-1) + ((j)-1) * ML]

static void ref_time_loop(double *cu, double *cv, double *z, double *h, double *u, double *unew,
                          double *uold, double *v, double *vnew, double *vold, double *p,
                          double *pnew, double *pold, int ML, int NL, int M, int N, double fsdx,
                          double fsdy, double dx, double dy, double alpha, double &tdt,
                          double &time, int itmax) {
  for (int it = 1; it <= itmax; ++it) {
    // compute_cu_cv_z_h
    for (int j = 1; j <= N; ++j)
      for (int i = 1; i <= M; ++i) {
        A(cu, i + 1, j) = 0.5 * (A(p, i + 1, j) + A(p, i, j)) * A(u, i + 1, j);
        A(cv, i, j + 1) = 0.5 * (A(p, i, j + 1) + A(p, i, j)) * A(v, i, j + 1);
        A(z, i + 1, j + 1) = (fsdx * (A(v, i + 1, j + 1) - A(v, i, j + 1)) -
                              fsdy * (A(u, i + 1, j + 1) - A(u, i + 1, j))) /
                             (A(p, i, j) + A(p, i + 1, j) + A(p, i + 1, j + 1) + A(p, i, j + 1));
        A(h, i, j) = A(p, i, j) + 0.25 * (A(u, i + 1, j) * A(u, i + 1, j) + A(u, i, j) * A(u, i, j) +
                                          A(v, i, j + 1) * A(v, i, j + 1) + A(v, i, j) * A(v, i, j));
      }
    // apply_bc_cu_cv_z_h
    for (int j = 1; j <= N; ++j) {
      A(cu, 1, j) = A(cu, ML, j);
      A(cv, ML, j + 1) = A(cv, 1, j + 1);
      A(z, 1, j + 1) = A(z, ML, j + 1);
      A(h, ML, j) = A(h, 1, j);
    }
    for (int i = 1; i <= M; ++i) {
      A(cu, i + 1, NL) = A(cu, i + 1, 1);
      A(cv, i, 1) = A(cv, i, NL);
      A(z, i + 1, 1) = A(z, i + 1, NL);
      A(h, i, NL) = A(h, i, 1);
    }
    A(cu, 1, NL) = A(cu, ML, 1);
    A(cv, ML, 1) = A(cv, 1, NL);
    A(z, 1, 1) = A(z, ML, NL);
    A(h, ML, NL) = A(h, 1, 1);

    const double tdts8 = tdt / 8.0, tdtsdx = tdt / dx, tdtsdy = tdt / dy;

    // update_unew_vnew_pnew
    for (int j = 1; j <= N; ++j)
      for (int i = 1; i <= M; ++i) {
        A(unew, i + 1, j) = A(uold, i + 1, j) +
                            tdts8 * (A(z, i + 1, j + 1) + A(z, i + 1, j)) *
                                (A(cv, i + 1, j + 1) + A(cv, i, j + 1) + A(cv, i, j) + A(cv, i + 1, j)) -
                            tdtsdx * (A(h, i + 1, j) - A(h, i, j));
        A(vnew, i, j + 1) = A(vold, i, j + 1) -
                            tdts8 * (A(z, i + 1, j + 1) + A(z, i, j + 1)) *
                                (A(cu, i + 1, j + 1) + A(cu, i, j + 1) + A(cu, i, j) + A(cu, i + 1, j)) -
                            tdtsdy * (A(h, i, j + 1) - A(h, i, j));
        A(pnew, i, j) = A(pold, i, j) - tdtsdx * (A(cu, i + 1, j) - A(cu, i, j)) -
                        tdtsdy * (A(cv, i, j + 1) - A(cv, i, j));
      }
    // apply_bc_unew_vnew_pnew
    for (int j = 1; j <= N; ++j) {
      A(unew, 1, j) = A(unew, ML, j);
      A(vnew, ML, j + 1) = A(vnew, 1, j + 1);
      A(pnew, ML, j) = A(pnew, 1, j);
    }
    for (int i = 1; i <= M; ++i) {
      A(unew, i + 1, NL) = A(unew, i + 1, 1);
      A(vnew, i, 1) = A(vnew, i, NL);
      A(pnew, i, NL) = A(pnew, i, 1);
    }
    A(unew, 1, NL) = A(unew, ML, 1);
    A(vnew, ML, 1) = A(vnew, 1, NL);
    A(pnew, ML, NL) = A(pnew, 1, 1);

    time += tdt;
    if (it > 1) {  // time_smoothing
      for (int j = 1; j <= NL; ++j)
        for (int i = 1; i <= ML; ++i) {
          A(uold, i, j) = A(u, i, j) + alpha * (A(unew, i, j) - 2.0 * A(u, i, j) + A(uold, i, j));
          A(vold, i, j) = A(v, i, j) + alpha * (A(vnew, i, j) - 2.0 * A(v, i, j) + A(vold, i, j));
          A(pold, i, j) = A(p, i, j) + alpha * (A(pnew, i, j) - 2.0 * A(p, i, j) + A(pold, i, j));
        }
    } else {  // first_cycle_copy
      tdt = tdt + tdt;
      for (int j = 1; j <= NL; ++j)
        for (int i = 1; i <= ML; ++i) {
          A(uold, i, j) = A(u, i, j);
          A(vold, i, j) = A(v, i, j);
          A(pold, i, j) = A(p, i, j);
        }
    }
    std::swap(u, unew);
    std::swap(v, vnew);
    std::swap(p, pnew);
  }
}

// Deterministic, non-periodic initial data (0-based i, j).
static double init_val(int k, int i, int j, int ML, int NL) {
  const double pi = 3.14159265358979323846;
  const double x = 2 * pi * i / ML, y = 2 * pi * j / NL;
  const double s = std::sin(12.9898 * i + 78.233 * j + 3.1 * k) * 43758.5453;
  const double noise = s - std::floor(s) - 0.5;  // in [-0.5, 0.5)
  switch (k) {
    case P: case PNEW: case POLD: return 5.0e4 + 2.0e3 * std::sin(x) * std::cos(y) + 100.0 * noise;
    case U: case UNEW: case UOLD: return 10.0 * std::sin(x + y) + 0.5 * noise;
    case V: case VNEW: case VOLD: return 8.0 * std::cos(x - y) + 0.5 * noise;
    default: return 777.0;  // overwritten everywhere if every point is covered
  }
}

static int run_case(int M, int N, int itmax) {
  const int ML = M + 1, NL = N + 1;
  const double dx = 1.0e5, dy = 1.0e5, alpha = 0.001;
  const double fsdx = 4.0 / dx, fsdy = 4.0 / dy;

  std::vector<std::vector<double>> ref(13, std::vector<double>(size_t(ML) * NL));
  DevView dev[13];
  HostView host[13];
  for (int k = 0; k < 13; ++k) {
    dev[k] = DevView(std::string(kNames[k]), ML, NL);
    host[k] = HostView(std::string("h_") + kNames[k], ML, NL);
    for (int j = 0; j < NL; ++j)
      for (int i = 0; i < ML; ++i) host[k](i, j) = ref[k][i + size_t(j) * ML] = init_val(k, i, j, ML, NL);
    Kokkos::deep_copy(dev[k], host[k]);
  }

  double tdt_ref = 90.0, time_ref = 0.0, tdt_k = 90.0, time_k = 0.0;
  ref_time_loop(ref[CU].data(), ref[CV].data(), ref[Z].data(), ref[H].data(), ref[U].data(),
                ref[UNEW].data(), ref[UOLD].data(), ref[V].data(), ref[VNEW].data(), ref[VOLD].data(),
                ref[P].data(), ref[PNEW].data(), ref[POLD].data(), ML, NL, M, N, fsdx, fsdy, dx, dy,
                alpha, tdt_ref, time_ref, itmax);

  kokkos_time_loop(dev[CU].data(), dev[CV].data(), dev[Z].data(), dev[H].data(), dev[U].data(),
                   dev[UNEW].data(), dev[UOLD].data(), dev[V].data(), dev[VNEW].data(),
                   dev[VOLD].data(), dev[P].data(), dev[PNEW].data(), dev[POLD].data(), ML, NL, M, N,
                   fsdx, fsdy, dx, dy, alpha, &tdt_k, &time_k, itmax);

  int bad = 0;
  double worst = 0.0;
  for (int k = 0; k < 13; ++k) {
    Kokkos::deep_copy(host[k], dev[k]);
    double maxabs = 0.0, maxdiff = 0.0;
    int wi = -1, wj = -1;
    for (int j = 0; j < NL; ++j)
      for (int i = 0; i < ML; ++i) {
        const double r = ref[k][i + size_t(j) * ML];
        maxabs = std::max(maxabs, std::fabs(r));
        const double d = std::fabs(host[k](i, j) - r);
        if (!(d <= maxdiff)) { maxdiff = d; wi = i; wj = j; }  // also catches NaN
      }
    const double rel = maxabs > 0 ? maxdiff / maxabs : maxdiff;
    worst = std::max(worst, rel);
    if (!(rel <= 1e-10)) {
      ++bad;
      std::printf("    MISMATCH %-5s max rel diff %.3e at (i=%d, j=%d) (0-based)\n", kNames[k], rel, wi, wj);
    }
  }
  const bool scalars_ok = (tdt_ref == tdt_k) && (time_ref == time_k);
  if (!scalars_ok) std::printf("    MISMATCH tdt/time: ref %g %g, kokkos %g %g\n", tdt_ref, time_ref, tdt_k, time_k);
  std::printf("  M=%d N=%d itmax=%d: %s (worst rel diff over 13 arrays %.2e)\n", M, N, itmax,
              (bad == 0 && scalars_ok) ? "PASS" : "FAIL", worst);
  return bad + (scalars_ok ? 0 : 1);
}

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  int fails = 0;
  {
    std::printf("Kokkos execution space: %s\n", SWM_EXEC_SPACE::name());
    fails += run_case(95, 47, 6);  // even: latest fields end up in u, v, p
    fails += run_case(95, 47, 7);  // odd:  latest fields end up in unew, vnew, pnew
    fails += run_case(5, 3, 4);    // tiny grid: edge threads coincide
    fails += run_case(64, 64, 1);  // first (forward) step only
  }
  Kokkos::finalize();
  std::printf(fails ? "FAILED\n" : "ALL PASSED\n");
  return fails ? 1 : 0;
}
