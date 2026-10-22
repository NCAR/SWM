// Timing driver for kokkos_time_loop.   usage: bench_kokkos M N iters   (iters even)
// Prints one RESULT line; the OpenACC driver (bench_acc.c) prints the same fields.

#include "shallow_kokkos.h"
#include "swm_init.h"

#include <Kokkos_Core.hpp>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>

using DevView = Kokkos::View<double **, Kokkos::LayoutLeft>;

int main(int argc, char **argv) {
  if (argc < 4) {
    std::fprintf(stderr, "usage: %s M N iters\n", argv[0]);
    return 2;
  }
  const int M = std::atoi(argv[1]), N = std::atoi(argv[2]);
  int iters = std::atoi(argv[3]);
  iters += iters % 2;  // even: the latest fields end up in u, v, p
  const int ML = M + 1, NL = N + 1;
  const double fsdx = 4.0 / SWM_DX, fsdy = 4.0 / SWM_DY;

  Kokkos::initialize(argc, argv);
  {
    const auto ti0 = std::chrono::steady_clock::now();
    DevView d[13];
    DevView::HostMirror h[13];
    for (int k = 0; k < 13; ++k) {
      d[k] = DevView(std::string("f") + std::to_string(k), ML, NL);
      h[k] = Kokkos::create_mirror_view(d[k]);
    }
    Kokkos::fence();
    const auto ti1 = std::chrono::steady_clock::now();  // device + host allocation done
    for (int k = 0; k < 13; ++k)
      for (int j = 0; j < NL; ++j)
        for (int i = 0; i < ML; ++i) h[k](i, j) = swm_init_val(k, i, j, M, N);
    const auto ti2 = std::chrono::steady_clock::now();  // host fill done
    double tdt = SWM_DT, t = 0.0;
    auto reset = [&]() {
      for (int k = 0; k < 13; ++k) Kokkos::deep_copy(d[k], h[k]);
      Kokkos::fence();
      tdt = SWM_DT;
      t = 0.0;
    };
    auto run = [&](int n) {
      kokkos_time_loop(d[SWM_CU].data(), d[SWM_CV].data(), d[SWM_Z].data(), d[SWM_H].data(),
                       d[SWM_U].data(), d[SWM_UNEW].data(), d[SWM_UOLD].data(), d[SWM_V].data(),
                       d[SWM_VNEW].data(), d[SWM_VOLD].data(), d[SWM_P].data(), d[SWM_PNEW].data(),
                       d[SWM_POLD].data(), ML, NL, M, N, fsdx, fsdy, SWM_DX, SWM_DY, SWM_ALPHA, &tdt, &t, n);
    };

    reset();  // host -> device copy of all 13 arrays
    const auto ti3 = std::chrono::steady_clock::now();
    {
      auto sec = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };
      std::printf("INIT M=%d N=%d alloc=%.4f host_fill=%.4f h2d_copy=%.4f total=%.4f\n", M, N,
                  sec(ti0, ti1), sec(ti1, ti2), sec(ti2, ti3), sec(ti0, ti3));
    }
    // warm-up disabled: time from a cold start, like the SWM programs and the Rust code
    // run(2);  // warm-up
    // reset();

    const auto t0 = std::chrono::steady_clock::now();
    run(iters);
    const auto t1 = std::chrono::steady_clock::now();
    const double sec = std::chrono::duration<double>(t1 - t0).count();

    // copy u, v, p back to the host (inside a second timer window: the SWM programs include this copy in their time)
    Kokkos::deep_copy(h[SWM_P], d[SWM_P]);
    Kokkos::deep_copy(h[SWM_U], d[SWM_U]);
    Kokkos::deep_copy(h[SWM_V], d[SWM_V]);
    const auto t2 = std::chrono::steady_clock::now();
    const double sec_copy = std::chrono::duration<double>(t2 - t0).count();

    double sum_p = 0.0, sum_u = 0.0, sum_v = 0.0;
    for (int j = 0; j < NL; ++j)
      for (int i = 0; i < ML; ++i) {
        sum_p += h[SWM_P](i, j);
        sum_u += h[SWM_U](i, j);
        sum_v += h[SWM_V](i, j);
      }
    // sec / ms_per_step: time loop only (as before); sec_copy / ms_per_step_copy: loop + copy of u, v, p to the host
    std::printf("RESULT impl=kokkos M=%d N=%d iters=%d sec=%.4f ms_per_step=%.4f Mpts_per_s=%.1f sum_p=%.12e sum_u=%.12e sec_copy=%.4f ms_per_step_copy=%.4f sum_v=%.12e\n",
                M, N, iters, sec, 1e3 * sec / iters, 1e-6 * double(M) * N * iters / sec, sum_p, sum_u,
                sec_copy, 1e3 * sec_copy / iters, sum_v);
  }
  Kokkos::finalize();
  return 0;
}
