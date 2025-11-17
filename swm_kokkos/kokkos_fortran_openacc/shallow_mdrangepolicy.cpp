/* Code converted from shallow_swap.c to use Kokkos for parallelism. */

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <Kokkos_Core.hpp>

#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))

#define M 256
#define N 256
#define M_LEN (M + 1)
#define N_LEN (N + 1)
#define ITMAX 4000
#define L_OUT true 
#define VAL_OUT false

using Layout = Kokkos::LayoutLeft;
using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemSpace = ExecSpace::memory_space;
using ViewMatrixType = Kokkos::View<double**, Layout, MemSpace, Kokkos::MemoryTraits<Kokkos::Restrict>>;
using HostViewMatrixType = Kokkos::View<double**, Layout, Kokkos::HostSpace>;

void write_to_file(auto array, int tM, int tN, const char *filename);
void print_to_file(auto array, int tM, int tN, const char *filename);

// C-style function prototype for the Fortran subroutine.
extern "C" {
    void fortran_time_loop(double* cu_raw_ptr, double* cv_raw_ptr, double* z_raw_ptr,
                          double* h_raw_ptr, double* u_raw_ptr, double* unew_raw_ptr,
                          double* uold_raw_ptr, double* v_raw_ptr, double* vnew_raw_ptr,
                          double* vold_raw_ptr, double* p_raw_ptr, double* pnew_raw_ptr,
                          double* pold_raw_ptr, int m_len, int n_len, int m, int n,
                          double fsdx, double fsdy, double dx, double dy, double alpha,
                          double& tdt, double& time, int itmax);
}

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
    double fsdx,fsdy;

    int mnmin;
  
    // timer variables
    double ctime,tcyc,time,ptime;

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
    Kokkos::parallel_for("init_psi_p", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}), KOKKOS_LAMBDA(const int i, const int j) {
      psi(i,j) = a * std::sin((i + .5) * di) * std::sin((j + .5) * dj);
      p(i,j) = pcf * (std::cos(2. * (i) * di) + std::cos(2. * (j) * dj)) + 50000.;
    });
    
    // Initialize velocities
    Kokkos::parallel_for("init_u_v", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M,N}), KOKKOS_LAMBDA(const int i, const int j) {
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

    Kokkos::parallel_for("init_old_arrays", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}), KOKKOS_LAMBDA(const int i, const int j) {
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

      mnmin = MIN(M,N);
      if constexpr (!std::is_same_v<MemSpace, Kokkos::HostSpace>) {
        Kokkos::deep_copy(u_host, u);
        Kokkos::deep_copy(v_host, v);
        Kokkos::deep_copy(p_host, p);
      }
      // printf(" initial diagonal elements of p\n");
      // for (int i=0; i<mnmin; i++) {
      //   printf("%f ",p_host(i,i));
      // }
      // printf("\n initial diagonal elements of u\n");
      // for (int i=0; i<mnmin; i++) {
      //   printf("%f ",u_host(i,i));
      // }
      // printf("\n initial diagonal elements of v\n");
      // for (int i=0; i<mnmin; i++) {
      //   printf("%f ",v_host(i,i));
      // }
      // printf("\n");
    }

    // Wait for the all the kernels above to complete before passing the raw pointers
    Kokkos::fence(); 

    // 5. Get the raw pointer from the Kokkos View
    double* cu_raw_ptr = cu.data();
    double* cv_raw_ptr = cv.data();
    double* z_raw_ptr = z.data();
    double* h_raw_ptr = h.data();
    double* u_raw_ptr = u.data();
    double* unew_raw_ptr = unew.data();
    double* uold_raw_ptr = uold.data();
    double* v_raw_ptr = v.data();
    double* vnew_raw_ptr = vnew.data();
    double* vold_raw_ptr = vold.data();
    double* p_raw_ptr = p.data();
    double* pnew_raw_ptr = pnew.data();
    double* pold_raw_ptr = pold.data();

    // Start timer
    Kokkos::Timer timer;
    time = 0.;

    // Call Fortran / OpenACC implementation for the time loop
    fortran_time_loop(cu_raw_ptr, cv_raw_ptr, z_raw_ptr, h_raw_ptr,
                      u_raw_ptr, unew_raw_ptr, uold_raw_ptr,
                      v_raw_ptr, vnew_raw_ptr, vold_raw_ptr,
                      p_raw_ptr, pnew_raw_ptr, pold_raw_ptr,
                      M_LEN, N_LEN, M, N, fsdx, fsdy,
                      dx, dy, alpha, tdt, time, ITMAX);

    // If ITMAX is odd, the valid data resides in the 'new' buffers
    // (unew, vnew, pnew). We must swap the C++ Views to match the
    // final state of the Fortran pointers.
    if (ITMAX % 2 == 1) {
      std::swap(u, unew);
      std::swap(v, vnew);
      std::swap(p, pnew);
    }

    // Try to use `if constexpr (!std::is_same_v<MemSpace, Kokkos::HostSpace>)` to use swap function
    //   for the host space, but somehow it is not compiled correctly for the device space.
    // Just use deep_copy for both spaces.
    Kokkos::deep_copy(u_host, u);
    Kokkos::deep_copy(v_host, v);
    Kokkos::deep_copy(p_host, p);

    // Output p, u, v fields and run times.q
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

      // printf(" diagonal elements of p\n");
      // for (int i=0; i<mnmin; i++) {
      //   printf("%f ",p_host(i,i));
      // }
      // printf("\n diagonal elements of u\n");
      // for (int i=0; i<mnmin; i++) {
      //   printf("%f ",u_host(i,i));
      // }
      // printf("\n diagonal elements of v\n");
      // for (int i=0; i<mnmin; i++) {
      //   printf("%f ",v_host(i,i));
      // }
      // printf("\n");

      ctime = timer.seconds();
      tcyc = ctime / ITMAX;
      printf(" cycle number %d total computer time %f time per cycle %f\n", ITMAX, ctime, tcyc);
    }
    /* 
      No need to free the raw pointers obtained from Kokkos::View::data().
      These pointers are managed by the Kokkos::View and do not require manual deallocation.
      Memory will be freed automatically when the View goes out of scope.
    */
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