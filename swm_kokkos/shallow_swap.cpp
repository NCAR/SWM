/* Code converted from shallow_swap.c to use Kokkos for parallelism. */

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <Kokkos_Core.hpp>

#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))

#define M 256
#define N 256
#define M_LEN (M + 1)
#define N_LEN (N + 1)
#define SIZE ((M_LEN)*(N_LEN))
#define ITMAX 4000
#define L_OUT true 
#define VAL_OUT false

using Layout = Kokkos::LayoutRight;
using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemSpace = ExecSpace::memory_space;
using ViewMatrixType = Kokkos::View<double**, Layout, MemSpace>;
using HostViewMatrixType = Kokkos::View<double**, Layout, Kokkos::HostSpace>;

void write_to_file(auto array, int tM, int tN, const char *filename);
void print_to_file(auto array, int tM, int tN, const char *filename);

int main(int argc, char **argv) {

  // Initialize Kokkos
  Kokkos::initialize( argc, argv );
  {

    // Allocate a large contiguous block for all variables
    constexpr int num_vars = 14;
    ViewMatrixType big_data("big_data", M_LEN * num_vars, N_LEN);

    // Create 2D views as subviews of big_data, and cast it to ViewMatrixType
    ViewMatrixType u    = Kokkos::subview(big_data, Kokkos::make_pair(0, M_LEN), Kokkos::make_pair(0, N_LEN));
    ViewMatrixType v    = Kokkos::subview(big_data, Kokkos::make_pair(M_LEN, 2 * M_LEN), Kokkos::make_pair(0, N_LEN));
    ViewMatrixType p    = Kokkos::subview(big_data, Kokkos::make_pair(2 * M_LEN, 3 * M_LEN), Kokkos::make_pair(0, N_LEN));
    ViewMatrixType unew = Kokkos::subview(big_data, Kokkos::make_pair(3 * M_LEN, 4 * M_LEN), Kokkos::make_pair(0, N_LEN));
    ViewMatrixType vnew = Kokkos::subview(big_data, Kokkos::make_pair(4 * M_LEN, 5 * M_LEN), Kokkos::make_pair(0, N_LEN));
    ViewMatrixType pnew = Kokkos::subview(big_data, Kokkos::make_pair(5 * M_LEN, 6 * M_LEN), Kokkos::make_pair(0, N_LEN));
    ViewMatrixType uold = Kokkos::subview(big_data, Kokkos::make_pair(6 * M_LEN, 7 * M_LEN), Kokkos::make_pair(0, N_LEN));
    ViewMatrixType vold = Kokkos::subview(big_data, Kokkos::make_pair(7 * M_LEN, 8 * M_LEN), Kokkos::make_pair(0, N_LEN));
    ViewMatrixType pold = Kokkos::subview(big_data, Kokkos::make_pair(8 * M_LEN, 9 * M_LEN), Kokkos::make_pair(0, N_LEN));
    ViewMatrixType cu   = Kokkos::subview(big_data, Kokkos::make_pair(9 * M_LEN, 10 * M_LEN), Kokkos::make_pair(0, N_LEN));
    ViewMatrixType cv   = Kokkos::subview(big_data, Kokkos::make_pair(10 * M_LEN, 11 * M_LEN), Kokkos::make_pair(0, N_LEN));
    ViewMatrixType z    = Kokkos::subview(big_data, Kokkos::make_pair(11 * M_LEN, 12 * M_LEN), Kokkos::make_pair(0, N_LEN));
    ViewMatrixType h    = Kokkos::subview(big_data, Kokkos::make_pair(12 * M_LEN, 13 * M_LEN), Kokkos::make_pair(0, N_LEN));
    ViewMatrixType psi  = Kokkos::subview(big_data, Kokkos::make_pair(13 * M_LEN, 14 * M_LEN), Kokkos::make_pair(0, N_LEN));

    double dt,tdt,dx,dy,a,alpha,el,pi;
    double tpi,di,dj,pcf;
    double tdts8,tdtsdx,tdtsdy,fsdx,fsdy;

    int mnmin,ncycle;
  
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

    // Start timer
    Kokkos::Timer timer;
    time = 0.;

    // ** Start of time loop ** 

    for (ncycle=1;ncycle<=ITMAX;++ncycle) {
      
      // Compute capital u, capital v, z and h

      // Compute cu
      Kokkos::parallel_for("compute_cu", Kokkos::RangePolicy<ExecSpace>(1, M_LEN), KOKKOS_LAMBDA(const int i) {
          for (int j = 0; j < N; ++j) {
            cu(i,j) = 0.5 * (p(i,j) + p(i-1,j)) * u(i,j);
          }
      });

      // Compute cv
      Kokkos::parallel_for("compute_cv", Kokkos::RangePolicy<ExecSpace>(0, M), KOKKOS_LAMBDA(const int i) {
          for (int j = 1; j < N_LEN; ++j) {
            cv(i,j) = 0.5 * (p(i,j) + p(i,j-1)) * v(i,j);
          }
      });

      // Compute z
      Kokkos::parallel_for("compute_z", Kokkos::RangePolicy<ExecSpace>(1, M_LEN), KOKKOS_LAMBDA(const int i) {
          for (int j = 1; j < N_LEN; ++j) {
            z(i,j) = (fsdx * (v(i,j) - v(i-1,j)) - fsdy * (u(i,j) - u(i,j-1))) / (p(i-1,j-1) + p(i,j-1) + p(i,j) + p(i-1,j));
          }
      });

      // Compute h
      Kokkos::parallel_for("compute_h", Kokkos::RangePolicy<ExecSpace>(0, M), KOKKOS_LAMBDA(const int i) {
          for (int j = 0; j < N; ++j) {
            h(i,j) = p(i,j) + 0.25 * ( u(i+1,j) * u(i+1,j) + u(i,j) * u(i,j) + v(i,j+1) * v(i,j+1) + v(i,j) * v(i,j) );
          }
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
      
      // Compute new values u,v and p

      tdts8 = tdt / 8.;
      tdtsdx = tdt / dx;
      tdtsdy = tdt / dy;

      // Compute unew
      Kokkos::parallel_for("compute_unew", Kokkos::RangePolicy<ExecSpace>(1, M_LEN), KOKKOS_LAMBDA(const int i) {
          for (int j = 0; j < N; ++j) {
            unew(i,j) = uold(i,j) + tdts8 * (z(i,j+1) + z(i,j)) * (cv(i,j+1) + cv(i-1,j+1) + cv(i-1,j) + cv(i,j)) - tdtsdx * (h(i,j) - h(i-1,j));
          }
      });

      // Compute vnew
      Kokkos::parallel_for("compute_vnew", Kokkos::RangePolicy<ExecSpace>(0, M), KOKKOS_LAMBDA(const int i) {
          for (int j = 1; j < N_LEN; ++j) {
            vnew(i,j) = vold(i,j) - tdts8 * (z(i+1,j) + z(i,j)) * (cu(i+1,j) + cu(i,j) + cu(i,j-1) + cu(i+1,j-1)) - tdtsdy * (h(i,j) - h(i,j-1));
          }
      });

      // Compute pnew
      Kokkos::parallel_for("compute_pnew", Kokkos::RangePolicy<ExecSpace>(0, M), KOKKOS_LAMBDA(const int i) {
          for (int j = 0; j < N; ++j) {
            pnew(i,j) = pold(i,j) - tdtsdx * (cu(i+1,j) - cu(i,j)) - tdtsdy * (cv(i,j+1) - cv(i,j));
          }
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

      time = time + dt;

      // Time smoothing and update for next cycle

      if ( ncycle > 1 ) {
        Kokkos::parallel_for("time_smoothing_u", Kokkos::RangePolicy<ExecSpace>(0,M_LEN), KOKKOS_LAMBDA(const int i) {
          for (int j = 0; j < N_LEN; ++j) {
            uold(i,j) = u(i,j) + alpha * (unew(i,j) - 2. * u(i,j) + uold(i,j));
          }
        });

        Kokkos::parallel_for("time_smoothing_v", Kokkos::RangePolicy<ExecSpace>(0,M_LEN), KOKKOS_LAMBDA(const int i) {
          for (int j = 0; j < N_LEN; ++j) {
            vold(i,j) = v(i,j) + alpha * (vnew(i,j) - 2. * v(i,j) + vold(i,j));
          }
        });

        Kokkos::parallel_for("time_smoothing_p", Kokkos::RangePolicy<ExecSpace>(0,M_LEN), KOKKOS_LAMBDA(const int i) {
          for (int j = 0; j < N_LEN; ++j) {
            pold(i,j) = p(i,j) + alpha * (pnew(i,j) - 2. * p(i,j) + pold(i,j));
          }
        }); 
      }
      else {
        tdt = tdt + tdt;
        Kokkos::parallel_for("first_cycle_copy", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}), KOKKOS_LAMBDA(const int i, const int j) {
            uold(i,j) = u(i,j);
            vold(i,j) = v(i,j);
            pold(i,j) = p(i,j);
        });
      }
      if constexpr (!std::is_same_v<MemSpace, Kokkos::HostSpace>) {
        Kokkos::fence();
        // Kokkos::parallel_for("swap_device_u_unew_views", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}), KOKKOS_LAMBDA(const int i, const int j) {
        //     auto temp = u(i,j);
        //     u(i,j) = unew(i,j);
        //     unew(i,j) = temp;
        // });
        // Kokkos::parallel_for("swap_device_v_vnew_views", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}), KOKKOS_LAMBDA(const int i, const int j) {
        //     auto temp = v(i,j);
        //     v(i,j) = vnew(i,j);
        //     vnew(i,j) = temp;
        // });
        // Kokkos::parallel_for("swap_device_p_pnew_views", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}), KOKKOS_LAMBDA(const int i, const int j) {
        //     auto temp = p(i,j);
        //     p(i,j) = pnew(i,j);
        //     pnew(i,j) = temp;
        // });
      }      
      // Swap the views
      std::swap(u, unew);
      std::swap(v, vnew);
      std::swap(p, pnew);
    } // ** End of time loop ** 

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