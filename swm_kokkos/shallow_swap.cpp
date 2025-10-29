/* Code converted from shallow_swap.c to use Kokkos for parallelism. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Kokkos_Core.hpp>

#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))

#define TRUE 1
#define FALSE 0
#define M 256
#define N 256
#define M_LEN (M + 1)
#define N_LEN (N + 1)
#define SIZE ((M_LEN)*(N_LEN))
#define ITMAX 4000
#define L_OUT TRUE
#define VAL_OUT TRUE 
//#define _COPY_

using ViewMatrixType = Kokkos::View<double**>;
void write_to_file(ViewMatrixType array, int tM, int tN, const char *filename);
void print_to_file(ViewMatrixType array, int tM, int tN, const char *filename);

int main(int argc, char **argv) {

  // Initialize Kokkos
  Kokkos::initialize( argc, argv );
  {

    // Define Kokkos Views
    ViewMatrixType u("u", M_LEN, N_LEN);
    ViewMatrixType v("v", M_LEN, N_LEN);
    ViewMatrixType p("p", M_LEN, N_LEN);
    ViewMatrixType unew("unew", M_LEN, N_LEN);
    ViewMatrixType vnew("vnew", M_LEN, N_LEN);
    ViewMatrixType pnew("pnew", M_LEN, N_LEN);
    ViewMatrixType uold("uold", M_LEN, N_LEN);
    ViewMatrixType vold("vold", M_LEN, N_LEN);
    ViewMatrixType pold("pold", M_LEN, N_LEN);
    ViewMatrixType cu("cu", M_LEN, N_LEN);
    ViewMatrixType cv("cv", M_LEN, N_LEN);
    ViewMatrixType z("z", M_LEN, N_LEN);
    ViewMatrixType h("h", M_LEN, N_LEN);
    ViewMatrixType psi("psi", M_LEN, N_LEN);

    double dt,tdt,dx,dy,a,alpha,el,pi;
    double tpi,di,dj,pcf;
    double tdts8,tdtsdx,tdtsdy,fsdx,fsdy;

    int mnmin,ncycle;
    int i,j;
  
    // timer variables 
    double mfs100,mfs200,mfs300;
    double t100,t200,t300;
    double tstart,ctime,tcyc,time,ptime;
    double c1,c2;

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
    Kokkos::parallel_for("init_psi_p", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}), KOKKOS_LAMBDA(const int i, const int j) {
        psi(i,j) = a * sin((i + .5) * di) * sin((j + .5) * dj);
        p(i,j) = pcf * (cos(2. * (i) * di) + cos(2. * (j) * dj)) + 50000.;
    });
    
    // Initialize velocities
    Kokkos::parallel_for("init_u_v", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {M,N}), KOKKOS_LAMBDA(const int i, const int j) {
        u(i+1,j) = -(psi(i+1,j+1) - psi(i+1,j)) / dy;
        v(i,j+1) = (psi(i+1,j+1) - psi(i,j+1)) / dx;
    });
      
    // Periodic continuation
    Kokkos::parallel_for("periodic_top_bottom_init", Kokkos::RangePolicy<>(0,N), KOKKOS_LAMBDA(const int j) {
      u(0,j) = u(M,j);
      v(M,j+1) = v(0,j+1);
    });

    Kokkos::parallel_for("periodic_left_right_init", Kokkos::RangePolicy<>(0,M), KOKKOS_LAMBDA(const int i) {
      u(i+1,N) = u(i+1,0);
      v(i,0) = v(i,N);
    });

    Kokkos::fence();

    u(0,N) = u(M,0);
    v(M,0) = v(0,N);
    
    Kokkos::parallel_for("init_old_arrays", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}), KOKKOS_LAMBDA(const int i, const int j) {
        unew(i,j) = u(i,j);
        vnew(i,j) = v(i,j);
        pnew(i,j) = p(i,j);
    });

    // Print initial values
    if ( L_OUT ) {
      Kokkos::fence();
      printf(" number of points in the x direction %d\n", N); 
      printf(" number of points in the y direction %d\n", M); 
      printf(" grid spacing in the x direction     %f\n", dx); 
      printf(" grid spacing in the y direction     %f\n", dy); 
      printf(" time step                           %f\n", dt); 
      printf(" time filter parameter               %f\n", alpha); 

      mnmin = MIN(M,N);
      printf(" initial diagonal elements of p\n");
      for (i=0; i<mnmin; i++) {
        printf("%f ",p(i,j));
      }
      printf("\n initial diagonal elements of u\n");
      for (i=0; i<mnmin; i++) {
        printf("%f ",u(i,j));
      }
      printf("\n initial diagonal elements of v\n");
      for (i=0; i<mnmin; i++) {
        printf("%f ",v(i,j));
      }
      printf("\n");
    }

    // Start timer
    Kokkos::Timer timer;

    // ** Start of time loop ** 

    for (ncycle=1;ncycle<=ITMAX;ncycle++) {
      
      // Compute capital u, capital v, z and h

      Kokkos::parallel_for("compute_cu_cv_z_h", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {M,N}), KOKKOS_LAMBDA(const int i, const int j) {
          cu(i+1,j) = .5 * (p(i+1,j) + p(i,j)) * u(i+1,j);
          cv(i,j+1) = .5 * (p(i,j+1) + p(i,j)) * v(i,j+1);
          z(i+1,j+1) = (fsdx * (v(i+1,j+1) - v(i,j+1)) - fsdy * (u(i+1,j+1) - u(i+1,j))) / (p(i,j) + p(i+1,j) + p(i+1,j+1) + p(i,j+1));
          h(i,j) = p(i,j) + .25 * (u(i+1,j) * u(i+1,j) + u(i,j) * u(i,j) + v(i,j+1) * v(i,j+1) + v(i,j) * v(i,j));
      });
  
      // Periodic continuation

      Kokkos::parallel_for("periodic_top_bottom_cu_cv_z_h", Kokkos::RangePolicy<>(0,N), KOKKOS_LAMBDA(const int j) {
        cu(0,j) = cu(M,j);
        cv(M,j+1) = cv(0,j+1);
        z(0,j+1) = z(M,j+1);
        h(M,j) = h(0,j);
      });

      Kokkos::parallel_for("periodic_left_right_cu_cv_z_h", Kokkos::RangePolicy<>(0,M), KOKKOS_LAMBDA(const int i) {
        cu(i+1,N) = cu(i+1,0);
        cv(i,0) = cv(i,N);
        z(i+1,0) = z(i+1,N);
        h(i,N) = h(i,0);
      });

      Kokkos::fence();
      cu(0,N) = cu(M,0);
      cv(M,0) = cv(0,N);
      z(0,0) = z(M,N);
      h(M,N) = h(0,0);
      
      // Compute new values u,v and p

      tdts8 = tdt / 8.;
      tdtsdx = tdt / dx;
      tdtsdy = tdt / dy;

      Kokkos::parallel_for("compute_unew_vnew_pnew", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {M,N}), KOKKOS_LAMBDA(const int i, const int j) {
          unew(i+1,j) = uold(i+1,j) + tdts8 * (z(i+1,j+1) + z(i+1,j)) * (cv(i+1,j+1) + cv(i,j+1) + cv(i,j) + cv(i+1,j)) - tdtsdx * (h(i+1,j) - h(i,j));
          vnew(i,j+1) = vold(i,j+1) - tdts8 * (z(i+1,j+1) + z(i,j+1)) * (cu(i+1,j+1) + cu(i,j+1) + cu(i,j) + cu(i+1,j)) - tdtsdy * (h(i,j+1) - h(i,j));
          pnew(i,j) = pold(i,j) - tdtsdx * (cu(i+1,j) - cu(i,j)) - tdtsdy * (cv(i,j+1) - cv(i,j));
      });

      // Periodic continuation

      Kokkos::parallel_for("periodic_top_bottom_unew_vnew_pnew", Kokkos::RangePolicy<>(0,N), KOKKOS_LAMBDA(const int j) {
        unew(0,j) = unew(M,j);
        vnew(M,j+1) = vnew(0,j+1);
        pnew(M,j) = pnew(0,j);
      });

      Kokkos::parallel_for("periodic_left_right_unew_vnew_pnew", Kokkos::RangePolicy<>(0,M), KOKKOS_LAMBDA(const int i) {
        unew(i+1,N) = unew(i+1,0);
        vnew(i,0) = vnew(i,N);
        pnew(i,N) = pnew(i,0);
      });

      Kokkos::fence();

      unew(0,N) = unew(M,0);
      vnew(M,0) = vnew(0,N);
      pnew(M,N) = pnew(0,0);

      time = time + dt;

      // Time smoothing and update for next cycle

      if ( ncycle > 1 ) {

        Kokkos::parallel_for("time_smoothing", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}), KOKKOS_LAMBDA(const int i, const int j) {
            uold(i,j) = u(i,j) + alpha * (unew(i,j) - 2. * u(i,j) + uold(i,j));
            vold(i,j) = v(i,j) + alpha * (vnew(i,j) - 2. * v(i,j) + vold(i,j));
            pold(i,j) = p(i,j) + alpha * (pnew(i,j) - 2. * p(i,j) + pold(i,j));
        });

        // Dependency
  #ifdef _COPY_
        Kokkos::parallel_for("copy_u_v_p_arrays", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}), KOKKOS_LAMBDA(const int i, const int j) {
            u(i,j) = unew(i,j);
            v(i,j) = vnew(i,j);
            p(i,j) = pnew(i,j);
        });
  #else
        // Swap the views
        Kokkos::fence();
        std::swap(u, unew);
        std::swap(v, vnew);
        std::swap(p, pnew);
  #endif     
      }
      else {
        tdt = tdt + tdt;

        Kokkos::parallel_for("first_cycle_copy", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {M_LEN,N_LEN}), KOKKOS_LAMBDA(const int i, const int j) {
            uold(i,j) = u(i,j);
            vold(i,j) = v(i,j);
            pold(i,j) = p(i,j);
        });

        Kokkos::fence();
        std::swap(unew, u);
        std::swap(vnew, v);
        std::swap(pnew, p);
      }
    } // ** End of time loop ** 

    // have to swap values back for printing

    Kokkos::fence();
    std::swap(u, unew);
    std::swap(v, vnew);
    std::swap(p, pnew);

    // Output p, u, v fields and run times.
    if(VAL_OUT){
      Kokkos::fence();
      write_to_file(pnew, M_LEN, N_LEN, "p.bin");
      write_to_file(unew, M_LEN, N_LEN, "u.bin");
      write_to_file(vnew, M_LEN, N_LEN, "v.bin");
      print_to_file(pnew, M_LEN, N_LEN, "p.txt");
      print_to_file(unew, M_LEN, N_LEN, "u.txt");
      print_to_file(vnew, M_LEN, N_LEN, "v.txt");
    }

    if (L_OUT) {
      Kokkos::fence();
      ptime = time / 3600.;

      printf(" cycle number %d model time in hours %f\n", ITMAX, ptime);
      printf(" diagonal elements of p\n");
      for (i=0; i<mnmin; i++) {
        printf("%f ",pnew(i,i));
      }
      printf("\n diagonal elements of u\n");
      for (i=0; i<mnmin; i++) {
        printf("%f ",unew(i,i));
      }
      printf("\n diagonal elements of v\n");
      for (i=0; i<mnmin; i++) {
        printf("%f ",vnew(i,i));
      }
      printf("\n");

      ctime = timer.seconds();
      tcyc = ctime / ITMAX;
      printf(" cycle number %d total computer time %f time per cycle %f\n", ITMAX, ctime, tcyc);
    }

  }
  
  // Finalize Kokkos
  Kokkos::finalize();

  return(0);
}

void print_to_file(ViewMatrixType array, int tM, int tN, const char *filename) {
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

void write_to_file(ViewMatrixType array, int tM, int tN, const char *filename) {
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
