/* Code converted from shallow_swap.c to use Kokkos for parallelism. */

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <Kokkos_Core.hpp>
#if defined(__has_include)
#if __has_include(<nvtx3/nvToolsExt.h>)
#include <nvtx3/nvToolsExt.h>
#define SWM_HAVE_NVTX 1
#endif
#endif
#ifndef SWM_HAVE_NVTX
#define nvtxRangePush(name) ((void)0)
#define nvtxRangePop() ((void)0)
#endif

#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))

#define M 2048
#define N M
#define M_LEN (M + 1)
#define N_LEN (N + 1)
#define SIZE ((M_LEN)*(N_LEN))
#define ITMAX 4000
#define L_OUT true
#define VAL_OUT false

using Layout = Kokkos::LayoutRight;
using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemSpace = ExecSpace::memory_space;
using ViewMatrixType = Kokkos::View<double**, Layout, MemSpace, Kokkos::MemoryTraits<Kokkos::Restrict>>;
using HostViewMatrixType = Kokkos::View<double**, Layout, Kokkos::HostSpace>;

void write_to_file(auto array, int tM, int tN, const char *filename);
void print_to_file(auto array, int tM, int tN, const char *filename);

int main(int argc, char **argv) {

  // Initialize Kokkos
  Kokkos::initialize( argc, argv );
  {
    // Number of variables to allocate
    constexpr int num_vars = 14;

    // Declare views before conditional so they are accessible later
    ViewMatrixType big_data;
    ViewMatrixType u, v, p, unew, vnew, pnew, uold, vold, pold, cu, cv, z, h, psi;

    if constexpr (std::is_same_v<Layout, Kokkos::LayoutLeft>) {
      // Allocate a large contiguous block for all variables
      big_data = ViewMatrixType("big_data", M_LEN, N_LEN * num_vars);
      // Create 2-D subviews for each variable
      u    = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(0 * N_LEN, 1 * N_LEN));
      v    = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(1 * N_LEN, 2 * N_LEN));
      p    = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(2 * N_LEN, 3 * N_LEN));
      unew = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(3 * N_LEN, 4 * N_LEN));
      vnew = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(4 * N_LEN, 5 * N_LEN));
      pnew = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(5 * N_LEN, 6 * N_LEN));
      uold = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(6 * N_LEN, 7 * N_LEN));
      vold = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(7 * N_LEN, 8 * N_LEN));
      pold = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(8 * N_LEN, 9 * N_LEN));
      cu   = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(9 * N_LEN, 10 * N_LEN));
      cv   = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(10 * N_LEN, 11 * N_LEN));
      z    = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(11 * N_LEN, 12 * N_LEN));
      h    = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(12 * N_LEN, 13 * N_LEN));
      psi  = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(13 * N_LEN, 14 * N_LEN));
    }
    else if constexpr (std::is_same_v<Layout, Kokkos::LayoutRight>) {
      big_data = ViewMatrixType("big_data", M_LEN * num_vars, N_LEN);
      u    = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(0, N_LEN));
      v    = Kokkos::subview(big_data, Kokkos::make_pair(M_LEN, 2 * M_LEN), Kokkos::make_pair(0, N_LEN));
      p    = Kokkos::subview(big_data, Kokkos::make_pair(2 * M_LEN, 3 * M_LEN), Kokkos::make_pair(0, N_LEN));
      unew = Kokkos::subview(big_data, Kokkos::make_pair(3 * M_LEN, 4 * M_LEN), Kokkos::make_pair(0, N_LEN));
      vnew = Kokkos::subview(big_data, Kokkos::make_pair(4 * M_LEN, 5 * M_LEN), Kokkos::make_pair(0, N_LEN));
      pnew = Kokkos::subview(big_data, Kokkos::make_pair(5 * M_LEN, 6 * M_LEN), Kokkos::make_pair(0, N_LEN));
      uold = Kokkos::subview(big_data, Kokkos::make_pair(6 * M_LEN, 7 * M_LEN), Kokkos::make_pair(0, N_LEN));
      vold = Kokkos::subview(big_data, Kokkos::make_pair(7 * M_LEN, 8 * M_LEN), Kokkos::make_pair(0, N_LEN));
      pold = Kokkos::subview(big_data, Kokkos::make_pair(8 * M_LEN, 9 * M_LEN), Kokkos::make_pair(0, N_LEN));
      cu   = Kokkos::subview(big_data, Kokkos::make_pair(9 * M_LEN, 10 * M_LEN), Kokkos::make_pair(0, N_LEN));
      cv   = Kokkos::subview(big_data, Kokkos::make_pair(10 * M_LEN, 11 * M_LEN), Kokkos::make_pair(0, N_LEN));
      z    = Kokkos::subview(big_data, Kokkos::make_pair(11 * M_LEN, 12 * M_LEN), Kokkos::make_pair(0, N_LEN));
      h    = Kokkos::subview(big_data, Kokkos::make_pair(12 * M_LEN, 13 * M_LEN), Kokkos::make_pair(0, N_LEN));
      psi  = Kokkos::subview(big_data, Kokkos::make_pair(13 * M_LEN, 14 * M_LEN), Kokkos::make_pair(0, N_LEN));
    }
    else {
      printf("Using unknown layout\n");
      return -1;
    }

    double dt,tdt,dx,dy,a,alpha,el,pi;
    double tpi,di,dj,pcf;
    double tdts8,tdtsdx,tdtsdy,fsdx,fsdy;

    int ncycle;
  
    // timer variables
    double ctime,tcyc,time,ptime;
    double t100 = 0., t200 = 0., t300 = 0.;

    // ** Initialisations ** 

    // Note below that two delta t (tdt) is set to dt on the first
    // cycle after which it is reset to dt+dt.
    dt = 90.;
    tdt = dt;
  
    dx = 100000.;
    dy = 100000.;
    fsdx = 4. / dx;
    fsdy = 4. / dy;

    a = 1000000.;
    alpha = .001;

    el = N * dx;
    pi = 4. * atan(1.);
    tpi = pi + pi;
    di = tpi / M;
    dj = tpi / N;
    pcf = pi * pi * a * a / (el * el);

    // Initial values of the stream function and p
    Kokkos::parallel_for("init_psi_p", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}, {1,256}), KOKKOS_LAMBDA(const int i, const int j) {
      psi(i,j) = a * std::sin((i + .5) * di) * std::sin((j + .5) * dj);
      p(i,j) = pcf * (std::cos(2. * (i) * di) + std::cos(2. * (j) * dj)) + 50000.;
    });
    
    // Initialize velocities
    Kokkos::parallel_for("init_u_v", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M,N}, {1,256}), KOKKOS_LAMBDA(const int i, const int j) {
      u(i+1,j) = -(psi(i+1,j+1) - psi(i+1,j)) / dy;
      v(i,j+1) = (psi(i+1,j+1) - psi(i,j+1)) / dx;
    });
      
    // Periodic continuation
    Kokkos::parallel_for("periodic_top_bottom_init", Kokkos::RangePolicy<ExecSpace>(0,N), KOKKOS_LAMBDA(const int j) {
      u(0,j) = u(M,j);
      v(M,j+1) = v(0,j+1);
    });

    Kokkos::parallel_for("periodic_left_right_init", Kokkos::RangePolicy<ExecSpace>(0,M), KOKKOS_LAMBDA(const int i) {
      u(i+1,N) = u(i+1,0);
      v(i,0) = v(i,N);
    });

    Kokkos::parallel_for("periodic_corners_init", Kokkos::RangePolicy<ExecSpace>(0,1), KOKKOS_LAMBDA(const int) {
      u(0,N) = u(M,0);
      v(M,0) = v(0,N);
    });

    Kokkos::parallel_for("init_old_arrays", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}, {1,256}), KOKKOS_LAMBDA(const int i, const int j) {
      uold(i,j) = u(i,j);
      vold(i,j) = v(i,j);
      pold(i,j) = p(i,j);
    });

    // Create host mirrors for output
    HostViewMatrixType u_host = Kokkos::create_mirror_view(u);
    HostViewMatrixType v_host = Kokkos::create_mirror_view(v);
    HostViewMatrixType p_host = Kokkos::create_mirror_view(p);

    // Print initial values
    if ( L_OUT ) {
      printf(" number of points in the x direction %d\n", N); 
      printf(" number of points in the y direction %d\n", M); 
      printf(" grid spacing in the x direction     %f\n", dx); 
      printf(" grid spacing in the y direction     %f\n", dy); 
      printf(" time step                           %f\n", dt); 
      printf(" time filter parameter               %f\n", alpha); 
    }

    // Start timer
    Kokkos::Timer timer;
    time = 0.;

    // ** Start of time loop ** 

    for (ncycle=1;ncycle<=ITMAX;++ncycle) {
      
      // Compute capital u, capital v, z and h
      nvtxRangePush("UpdateIntermediateVariables");
      Kokkos::Timer timer100;

      // Compute cu, cv, z, h (fused)
      Kokkos::parallel_for("compute_cu_cv_z_h", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M,N}, {1,256}), KOKKOS_LAMBDA(const int i, const int j) {
        double p_ij = p(i,j);
        double p_i1j = p(i+1,j);
        double p_ij1 = p(i,j+1);
        double p_i1j1 = p(i+1,j+1);
        double u_ij = u(i,j);
        double u_i1j = u(i+1,j);
        double u_i1j1 = u(i+1,j+1);
        double v_ij = v(i,j);
        double v_ij1 = v(i,j+1);
        double v_i1j1 = v(i+1,j+1);
        cu(i+1,j) = 0.5 * (p_i1j + p_ij) * u_i1j;
        cv(i,j+1) = 0.5 * (p_ij1 + p_ij) * v_ij1;
        z(i+1,j+1) = (fsdx * (v_i1j1 - v_ij1) - fsdy * (u_i1j1 - u_i1j)) / (p_ij + p_i1j + p_i1j1 + p_ij1);
        h(i,j) = p_ij + 0.25 * (u_i1j * u_i1j + u_ij * u_ij + v_ij1 * v_ij1 + v_ij * v_ij);
      });
  
      // Periodic continuation
      Kokkos::parallel_for("periodic_top_bottom_cu_cv_z_h", Kokkos::RangePolicy<ExecSpace>(0,N), KOKKOS_LAMBDA(const int j) {
        cu(0,j) = cu(M,j);
        cv(M,j+1) = cv(0,j+1);
        z(0,j+1) = z(M,j+1);
        h(M,j) = h(0,j);
      });

      Kokkos::parallel_for("periodic_left_right_cu_cv_z_h", Kokkos::RangePolicy<ExecSpace>(0,M), KOKKOS_LAMBDA(const int i) {
        cu(i+1,N) = cu(i+1,0);
        cv(i,0) = cv(i,N);
        z(i+1,0) = z(i+1,N);
        h(i,N) = h(i,0);
      });

      Kokkos::parallel_for("periodic_corner_cu_cv_z_h", Kokkos::RangePolicy<ExecSpace>(0, 1), KOKKOS_LAMBDA(const int) {
        cu(0,N) = cu(M,0);
        cv(M,0) = cv(0,N);
        z(0,0) = z(M,N);
        h(M,N) = h(0,0);
      });
      
      Kokkos::fence();
      t100 += timer100.seconds();
      nvtxRangePop();

      // Compute new values u,v and p
      nvtxRangePush("UpdateNewVariables");
      Kokkos::Timer timer200;

      tdts8 = tdt / 8.;
      tdtsdx = tdt / dx;
      tdtsdy = tdt / dy;

      // Compute unew, vnew, pnew (fused)
      Kokkos::parallel_for("compute_unew_vnew_pnew", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M,N}, {1,256}), KOKKOS_LAMBDA(const int i, const int j) {
        double z_i1j1 = z(i+1,j+1);
        double z_i1j = z(i+1,j);
        double z_ij1 = z(i,j+1);
        double cv_i1j1 = cv(i+1,j+1);
        double cv_ij1 = cv(i,j+1);
        double cv_ij = cv(i,j);
        double cv_i1j = cv(i+1,j);
        double cu_i1j1 = cu(i+1,j+1);
        double cu_ij1 = cu(i,j+1);
        double cu_ij = cu(i,j);
        double cu_i1j = cu(i+1,j);
        double h_i1j = h(i+1,j);
        double h_ij = h(i,j);
        double h_ij1 = h(i,j+1);
        unew(i+1,j) = uold(i+1,j) + tdts8 * (z_i1j1 + z_i1j) * (cv_i1j1 + cv_ij1 + cv_ij + cv_i1j) - tdtsdx * (h_i1j - h_ij);
        vnew(i,j+1) = vold(i,j+1) - tdts8 * (z_i1j1 + z_ij1) * (cu_i1j1 + cu_ij1 + cu_ij + cu_i1j) - tdtsdy * (h_ij1 - h_ij);
        pnew(i,j) = pold(i,j) - tdtsdx * (cu_i1j - cu_ij) - tdtsdy * (cv_ij1 - cv_ij);
      });

      // Periodic continuation
      Kokkos::parallel_for("periodic_top_bottom_unew_vnew_pnew", Kokkos::RangePolicy<ExecSpace>(0,N), KOKKOS_LAMBDA(const int j) {
        unew(0,j) = unew(M,j);
        vnew(M,j+1) = vnew(0,j+1);
        pnew(M,j) = pnew(0,j);
      });

      Kokkos::parallel_for("periodic_left_right_unew_vnew_pnew", Kokkos::RangePolicy<ExecSpace>(0,M), KOKKOS_LAMBDA(const int i) {
        unew(i+1,N) = unew(i+1,0);
        vnew(i,0) = vnew(i,N);
        pnew(i,N) = pnew(i,0);
      });

      Kokkos::parallel_for("periodic_corner_unew_vnew_pnew", Kokkos::RangePolicy<ExecSpace>(0, 1), KOKKOS_LAMBDA(const int) {
        unew(0,N) = unew(M,0);
        vnew(M,0) = vnew(0,N);
        pnew(M,N) = pnew(0,0);
      });

      Kokkos::fence();
      t200 += timer200.seconds();
      nvtxRangePop();

      time = time + dt;

      // Time smoothing and update for next cycle
      nvtxRangePush("UpdateOldVariables");
      Kokkos::Timer timer300;

      if ( ncycle > 1 ) {
        Kokkos::parallel_for("time_smoothing_uvp", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}, {1,256}), KOKKOS_LAMBDA(const int i, const int j) {
          double u_ij = u(i,j);
          double v_ij = v(i,j);
          double p_ij = p(i,j);
          double uold_ij = uold(i,j);
          double vold_ij = vold(i,j);
          double pold_ij = pold(i,j);
          uold(i,j) = u_ij + alpha * (unew(i,j) - 2. * u_ij + uold_ij);
          vold(i,j) = v_ij + alpha * (vnew(i,j) - 2. * v_ij + vold_ij);
          pold(i,j) = p_ij + alpha * (pnew(i,j) - 2. * p_ij + pold_ij);
        });
      }
      else {
        tdt = tdt + tdt;
        Kokkos::parallel_for("first_cycle_copy", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}, {1,256}), KOKKOS_LAMBDA(const int i, const int j) {
          uold(i,j) = u(i,j);
          vold(i,j) = v(i,j);
          pold(i,j) = p(i,j);
        });
      }

        Kokkos::fence();
      t300 += timer300.seconds();
      nvtxRangePop();

      // Swap the views
      std::swap(u, unew);
      std::swap(v, vnew);
      std::swap(p, pnew);
    } // ** End of time loop ** 

    ctime = timer.seconds();

    // Try to use `if constexpr (!std::is_same_v<MemSpace, Kokkos::HostSpace>)` to use swap function
    //   for the host space, but somehow it is not compiled correctly for the device space.
    // Just use deep_copy for both spaces.
    Kokkos::deep_copy(u_host, u);
    Kokkos::deep_copy(v_host, v);
    Kokkos::deep_copy(p_host, p);

    // Output p, u, v fields and run times.
    if(VAL_OUT) {
      write_to_file(p_host, M_LEN, N_LEN, "p.bin");
      write_to_file(u_host, M_LEN, N_LEN, "u.bin");
      write_to_file(v_host, M_LEN, N_LEN, "v.bin");
      print_to_file(p_host, M_LEN, N_LEN, "p.txt");
      print_to_file(u_host, M_LEN, N_LEN, "u.txt");
      print_to_file(v_host, M_LEN, N_LEN, "v.txt");
    }

    if (L_OUT) {
      ptime = time / 3600.;
      printf(" cycle number %d model time in hours %f\n", ITMAX, ptime);

      tcyc = ctime / ITMAX;
      printf(" cycle number %d total computer time %f time per cycle %f\n", ITMAX, ctime, tcyc);

      double mfs100 = double(ITMAX) * double(M) * double(N) * 16.e-6 / t100;
      double mfs200 = double(ITMAX) * double(M) * double(N) * 26.e-6 / t200;
      double mfs300 = double(ITMAX-1) * double(M_LEN) * double(N_LEN) * 15.e-6 / t300;
      printf(" time and megaflops for loop 100 %f %f\n", t100, mfs100);
      printf(" time and megaflops for loop 200 %f %f\n", t200, mfs200);
      printf(" time and megaflops for loop 300 %f %f\n", t300, mfs300);
    }

  }
  // Finalize Kokkos
  Kokkos::finalize();

  return(0);
}

void print_to_file(auto array, int tM, int tN, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        return;
    }
    for (int i = 0; i < tM; i++) {
        for (int j = 0; j < tN; j++) {
            fprintf(file, "%f ", array(i,j));
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void write_to_file(auto array, int tM, int tN, const char *filename) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        return;
    }
    for (int i = 0; i < tM; i++) {
        for (int j = 0; j < tN; j++) {
          fwrite(&array(i,j), sizeof(double), 1, file);
        }
    }
    fclose(file);
}