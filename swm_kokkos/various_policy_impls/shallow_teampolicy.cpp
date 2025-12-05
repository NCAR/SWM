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
using teamPolicy = Kokkos::TeamPolicy<ExecSpace>;
using memberType = teamPolicy::member_type;

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
    Kokkos::parallel_for("init_psi_p", teamPolicy(M_LEN, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
        const int i = team.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N_LEN), [&](const int j) {
            psi(i,j) = a * std::sin((i + .5) * di) * std::sin((j + .5) * dj);
            p(i,j) = pcf * (std::cos(2. * (i) * di) + std::cos(2. * (j) * dj)) + 50000.;
        });
    });
    
    // Initialize velocities
    Kokkos::parallel_for("init_u_v", teamPolicy(M, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
        const int i = team.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int j) {
            u(i+1,j) = -(psi(i+1,j+1) - psi(i+1,j)) / dy;
            v(i,j+1) = (psi(i+1,j+1) - psi(i,j+1)) / dx;
        });
    });
      
    // Periodic continuation
    Kokkos::parallel_for("periodic_top_bottom_init", teamPolicy(N, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
      const int j = team.league_rank();
      u(0,j) = u(M,j);
      v(M,j+1) = v(0,j+1);
    });

    Kokkos::parallel_for("periodic_left_right_init", teamPolicy(M, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
      const int i = team.league_rank();
      u(i+1,N) = u(i+1,0);
      v(i,0) = v(i,N);
    });

    Kokkos::parallel_for("periodic_corners_init", teamPolicy(1, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
      u(0,N) = u(M,0);
      v(M,0) = v(0,N);
    });

    Kokkos::parallel_for("init_old_arrays", teamPolicy(M_LEN, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N_LEN), [&](const int j) {
        uold(i,j) = u(i,j);
        vold(i,j) = v(i,j);
        pold(i,j) = p(i,j);
      });
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
      Kokkos::parallel_for("compute_cu", teamPolicy(M, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
          const int i = team.league_rank() + 1;
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int j) {
            cu(i,j) = 0.5 * (p(i,j) + p(i-1,j)) * u(i,j);
          });
      });

      // Compute cv
      Kokkos::parallel_for("compute_cv", teamPolicy(M, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
          const int i = team.league_rank();
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 1, N_LEN), [&](const int j) {
            cv(i,j) = 0.5 * (p(i,j) + p(i,j-1)) * v(i,j);
          });
      });

      // Compute z
      Kokkos::parallel_for("compute_z", teamPolicy(M, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
          const int i = team.league_rank() + 1;
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 1, N_LEN), [&](const int j) {
            z(i,j) = (fsdx * (v(i,j) - v(i-1,j)) - fsdy * (u(i,j) - u(i,j-1))) / (p(i-1,j-1) + p(i,j-1) + p(i,j) + p(i-1,j));
          });
      });

      // Compute h
      Kokkos::parallel_for("compute_h", teamPolicy(M, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
          const int i = team.league_rank();
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int j) {
            h(i,j) = p(i,j) + 0.25 * (u(i+1,j) * u(i+1,j) + u(i,j) * u(i,j) + v(i,j+1) * v(i,j+1) + v(i,j) * v(i,j));
          });
      });
  
      // Periodic continuation
      Kokkos::parallel_for("periodic_top_bottom_cu_cv_z_h", teamPolicy(N, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
        const int j = team.league_rank();
        cu(0,j) = cu(M,j);
        cv(M,j+1) = cv(0,j+1);
        z(0,j+1) = z(M,j+1);
        h(M,j) = h(0,j);
      });

      Kokkos::parallel_for("periodic_left_right_cu_cv_z_h", teamPolicy(M, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
        const int i = team.league_rank();
        cu(i+1,N) = cu(i+1,0);
        cv(i,0) = cv(i,N);
        z(i+1,0) = z(i+1,N);
        h(i,N) = h(i,0);
      });

      Kokkos::parallel_for("periodic_corner_cu_cv_z_h", teamPolicy(1, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
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
      Kokkos::parallel_for("compute_unew", teamPolicy(M, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
          const int i = team.league_rank() + 1;
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int j) {
            unew(i,j) = uold(i,j) + tdts8 * (z(i,j+1) + z(i,j)) * (cv(i,j+1) + cv(i-1,j+1) + cv(i-1,j) + cv(i,j)) - tdtsdx * (h(i,j) - h(i-1,j));
          });
      });

      // Compute vnew
      Kokkos::parallel_for("compute_vnew", teamPolicy(M, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
          const int i = team.league_rank();
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 1, N_LEN), [&](const int j) {
            vnew(i,j) = vold(i,j) - tdts8 * (z(i+1,j) + z(i,j)) * (cu(i+1,j) + cu(i,j) + cu(i,j-1) + cu(i+1,j-1)) - tdtsdy * (h(i,j) - h(i,j-1));
          });
      });

      // Compute pnew
      Kokkos::parallel_for("compute_pnew", teamPolicy(M, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
          const int i = team.league_rank();
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int j) {
            pnew(i,j) = pold(i,j) - tdtsdx * (cu(i+1,j) - cu(i,j)) - tdtsdy * (cv(i,j+1) - cv(i,j));
          });
      });

      // Periodic continuation
      Kokkos::parallel_for("periodic_top_bottom_unew_vnew_pnew", teamPolicy(N, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
        const int j = team.league_rank();
        unew(0,j) = unew(M,j);
        vnew(M,j+1) = vnew(0,j+1);
        pnew(M,j) = pnew(0,j);
      });

      Kokkos::parallel_for("periodic_left_right_unew_vnew_pnew", teamPolicy(M, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
        const int i = team.league_rank();
        unew(i+1,N) = unew(i+1,0);
        vnew(i,0) = vnew(i,N);
        pnew(i,N) = pnew(i,0);
      });

      Kokkos::parallel_for("periodic_corner_unew_vnew_pnew", teamPolicy(1, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
        unew(0,N) = unew(M,0);
        vnew(M,0) = vnew(0,N);
        pnew(M,N) = pnew(0,0);
      });

      time = time + dt;

      // Time smoothing and update for next cycle

      if ( ncycle > 1 ) {
        Kokkos::parallel_for("time_smoothing_u", teamPolicy(M_LEN, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
          const int i = team.league_rank();
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N_LEN), [&](const int j) {
            uold(i,j) = u(i,j) + alpha * (unew(i,j) - 2. * u(i,j) + uold(i,j));
          });
        });

        Kokkos::parallel_for("time_smoothing_v", teamPolicy(M_LEN, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
          const int i = team.league_rank();
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N_LEN), [&](const int j) {
            vold(i,j) = v(i,j) + alpha * (vnew(i,j) - 2. * v(i,j) + vold(i,j));
          });
        });

        Kokkos::parallel_for("time_smoothing_p", teamPolicy(M_LEN, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
          const int i = team.league_rank();
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N_LEN), [&](const int j) {
            pold(i,j) = p(i,j) + alpha * (pnew(i,j) - 2. * p(i,j) + pold(i,j));
          });
        });
      }
      else {
        tdt = tdt + tdt;
        Kokkos::parallel_for("first_cycle_copy", teamPolicy(M_LEN, Kokkos::AUTO), KOKKOS_LAMBDA(const memberType& team) {
            const int i = team.league_rank();
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N_LEN), [&](const int j) {
              uold(i,j) = u(i,j);
              vold(i,j) = v(i,j);
              pold(i,j) = p(i,j);
            });
        });
      }
      // My test shows it is better to use fence and then swap for device views, rather than doing a device copy without fence.
      if constexpr (!std::is_same_v<MemSpace, Kokkos::HostSpace>) {
        Kokkos::fence();
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