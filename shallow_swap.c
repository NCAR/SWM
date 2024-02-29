//#define _COPY_

/* Code converted from shallow_base.f90 using F2C-ACC program. 
 * Manually replaced: 
 * - WRITE statements with printf
 * - MOD operator with % 
 * - system_clock with wtime
 * Fixed several of the array references which had x dimension as 1, 
 * instead of M_LEN. 
 * Fixed values set using d and e notation. 
 * (7 June 2011)
 ***************
 * 'Pure' C version developed by G.D Riley (UoM) (25 Jan 2012)
 * removed all ftocmacros
 * used sin and cos not sinf and cosf (since all data are doubles)
 * needed to declare arrays +1 to cope with Fortran indexing
 * Compile:
 * gcc -O2 -c wtime.c
 * gcc -O2 -o sb shallow_base.c -lm wtime.o
 * May need to set 'ulimit -s unlimited' to run large problems (e.g. 512x512)
 * Results are consistent with Fortran version of the code
 *
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENACC
#include <openacc.h>
#endif
#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))

#define TRUE 1
#define FALSE 0
#define M 16
#define N 16
#define M_LEN (M + 1)
#define N_LEN (N + 1)
#define SIZE ((M_LEN)*(N_LEN))
#define ITMAX 1
#define L_OUT TRUE
#define VAL_OUT TRUE

extern double wtime(); 
extern void dswap(double **a, double **b);
extern void write_to_file(double *array, int tM, int tN, const char *filename);
extern void print_to_file(double *array, int tM, int tN, const char *filename);

//! Benchmark weather prediction program for comparing the
//! preformance of current supercomputers. The model is
//! based on the paper - The Dynamics of Finite-Difference
//! Models of the Shallow-Water Equations, by Robert Sadourny
//! J. Atm. Sciences, Vol 32, No 4, April 1975.
//!     
//! Code by Paul N. Swarztrauber, National Center for
//! Atmospheric Research, Boulder, Co,  October 1984.
//! Modified by Juliana Rew, NCAR, January 2006
//!
//! In this version, shallow4.f, initial and calculated values
//! of U, V, and P are written to a netCDF file
//! for later use in visualizing the results. The netCDF data
//! management library is freely available from
//! http://www.unidata.ucar.edu/software/netcdf
//! This code is still serial but has been brought up to modern
//! Fortran constructs and uses portable intrinsic Fortran 90 timing routines. 
//! This can be compiled on the IBM SP using:
//! xlf90 -qmaxmem=-1 -g -o shallow4 -qfixed=132 -qsclk=micro \
//! -I/usr/local/include shallow4.f -L/usr/local/lib32/r4i4 -l netcdf
//! where the -L and -I point to local installation of netCDF
//!     
//! Changes from shallow4.f (Annette Osprey, January 2010):
//! - Converted to free-form fortran 90.  
//! - Some tidying up of old commented-out code.   
//! - Explicit type declarations.
//! - Variables n, m, ITMAX and mprint read in from namelist. 
//! - Dynamic array allocation.
//! - Only write to netcdf at mprint timesteps.
//! - Don't write wrap-around points to NetCDF file.
//! - Use 8-byte reals.
//!
//! Further changes (Annette Osprey & Graham Riley, February 2011): 
//! - Remove unnecessary halo updates.
//! - Split loops to improve TLB access.
//! - Modify timers to accumulate loop times over whole run. 
//! - Remove old-style indentation. 
//!
//! Minimal serial version (26 May 2011)

int main(int argc, char **argv) {
  
  // solution arrays
  double *u,*v,*p,*unew,*vnew,*pnew,*uold,*vold,*pold,*cu,*cv,*z,*h,*psi;

  u = (double *)malloc(sizeof(double)*M_LEN*N_LEN);
  v = (double *)malloc(sizeof(double)*M_LEN*N_LEN);
  p = (double *)malloc(sizeof(double)*M_LEN*N_LEN);
  unew = (double *)malloc(sizeof(double)*M_LEN*N_LEN);
  vnew = (double *)malloc(sizeof(double)*M_LEN*N_LEN);
  pnew = (double *)malloc(sizeof(double)*M_LEN*N_LEN);
  uold = (double *)malloc(sizeof(double)*M_LEN*N_LEN);
  vold = (double *)malloc(sizeof(double)*M_LEN*N_LEN);
  pold = (double *)malloc(sizeof(double)*M_LEN*N_LEN);
  cu = (double *)malloc(sizeof(double)*M_LEN*N_LEN);
  cv = (double *)malloc(sizeof(double)*M_LEN*N_LEN);
  z = (double *)malloc(sizeof(double)*M_LEN*N_LEN);
  h = (double *)malloc(sizeof(double)*M_LEN*N_LEN);
  psi = (double *)malloc(sizeof(double)*M_LEN*N_LEN);

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

  for (i=0;i<M_LEN;i++) {
    for (j=0;j<N_LEN;j++) {
      int idx = i*N_LEN+j;
      psi[idx] = a * sin((i + .5) * di) * sin((j + .5) * dj);
      p[idx] = pcf * (cos(2. * (i) * di) + cos(2. * (j) * dj)) + 50000.;
    }
  }
  
  // Initialize velocities

  for (i=0;i<M;i++) {
    for (j=0;j<N;j++) {
      int idx01 = (i*N_LEN) + j+1; //[i][j+1]
      int idx10 = ((i+1)*N_LEN) + j; //[i+1][j]
      int idx11 = ((i+1)*N_LEN) + j+1; //[i+1][j+1]
      u[idx10] = -(psi[idx11] - psi[idx10]) / dy;
      v[idx01] = (psi[idx11] - psi[idx01]) / dx;
    }
  }
     
  // Periodic continuation

  for (j=0;j<N;j++) {
    u[j] = u[M*N_LEN+j]; 
    v[M*N_LEN+j+1] = v[j + 1];
  }

  for (i=0;i<M;i++) {
    u[(i + 1)*N_LEN +N] = u[(i + 1)*N_LEN];
    v[i*N_LEN] = v[(i*N_LEN)+N];
  }

  u[N] = u[M*N_LEN];
  v[M*N_LEN] = v[N];

  for (i=0;i<M_LEN;i++) {
    for (j=0;j<N_LEN;j++) {
      int idx = i*N_LEN+j;
      uold[idx] = u[idx];
      vold[idx] = v[idx];
      pold[idx] = p[idx];
    }
  }

  // Print initial values
  if ( L_OUT ) {
    printf(" number of points in the x direction %d\n", N); 
    printf(" number of points in the y direction %d\n", M); 
    printf(" grid spacing in the x direction     %f\n", dx); 
    printf(" grid spacing in the y direction     %f\n", dy); 
    printf(" time step                           %f\n", dt); 
    printf(" time filter parameter               %f\n", alpha); 

    mnmin = MIN(M,N);
    printf(" initial diagonal elements of p\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",p[i*N_LEN+i]);
    }
    printf("\n initial diagonal elements of u\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",u[i*N_LEN+i]);
    }
    printf("\n initial diagonal elements of v\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",v[i*N_LEN+i]);
    }
    printf("\n");
  }

  // Start timer
  tstart = wtime(); 
  time = 0.;
  t100 = 0.;
  t200 = 0.;
  t300 = 0.;

  // ** Start of time loop ** 

  for (ncycle=1;ncycle<=ITMAX;ncycle++) {
    
    // Compute capital u, capital v, z and h
    c1 = wtime();  

    for (i=0;i<M;i++) {
      for (j=0;j<N;j++) {
        int idx00 = (i*N_LEN) + j;
        int idx01 = (i*N_LEN) + j+1;
        int idx10 = ((i+1)*N_LEN) + j;
        int idx11 = ((i+1)*N_LEN) + j+1;
        cu[idx10] = .5 * (p[idx10] + p[idx00]) * u[idx10];
        cv[idx01] = .5 * (p[idx01] + p[idx00]) * v[idx01];
        z[idx11] = (fsdx * (v[idx11] - v[idx01]) - fsdy * (u[idx11] - u[idx10])) / (p[idx00] + p[idx10] + p[idx11] + p[idx01]);
        h[idx00] = p[idx00] + .25 * (u[idx10] * u[idx10] + u[idx00] * u[idx00] + v[idx01] * v[idx01] + v[idx00] * v[idx00]);
      }
    }

    c2 = wtime();  
    t100 = t100 + (c2 - c1); 

    // Periodic continuation

    for (j=0;j<N;j++) {
      cu[j] = cu[M*N_LEN + j];
      cv[M*N_LEN + j + 1] = cv[j + 1];
      z[j + 1] = z[M*N_LEN + j + 1];
      h[M*N_LEN +j] = h[j];
    }

    for (i=0;i<M;i++) {
      cu[(i + 1)*N_LEN +N] = cu[(i + 1)*N_LEN];
      cv[i*N_LEN] = cv[i*N_LEN + N];
      z[(i + 1)*N_LEN] = z[(i+1)*N_LEN + N];
      h[i*N_LEN + N] = h[i*N_LEN];
    }

    cu[N] = cu[M*N_LEN];
    cv[M*N_LEN] = cv[N];
    z[0] = z[M*N_LEN+N];
    h[M*N_LEN+N] = h[0];
     
    // Compute new values u,v and p

    tdts8 = tdt / 8.;
    tdtsdx = tdt / dx;
    tdtsdy = tdt / dy;

    c1 = wtime(); 

    for (i=0;i<M;i++) {
      for (j=0;j<N;j++) {
        int idx00 = (i*N_LEN) + j;
        int idx01 = (i*N_LEN) + j+1;
        int idx10 = ((i+1)*N_LEN) + j;
        int idx11 = ((i+1)*N_LEN) + j+1;
        unew[idx10] = uold[idx10] + tdts8 * (z[idx11] + z[idx10]) * (cv[idx11] + cv[idx01] + cv[idx00] + cv[idx10]) - tdtsdx * (h[idx10] - h[idx00]);
        vnew[idx01] = vold[idx01] - tdts8 * (z[idx11] + z[idx01]) * (cu[idx11] + cu[idx01] + cu[idx00] + cu[idx10]) - tdtsdy * (h[idx01] - h[idx00]);
        pnew[idx00] = pold[idx00] - tdtsdx * (cu[idx10] - cu[idx00]) - tdtsdy * (cv[idx01] - cv[idx00]);
      }
    }

    c2 = wtime();  
    t200 = t200 + (c2 - c1); 

    // Periodic continuation

    for (j=0;j<N;j++) {
      unew[j] = unew[M*N_LEN+j];
      vnew[M*N_LEN +j + 1] = vnew[j + 1];
      pnew[M*N_LEN +j] = pnew[j];
    }

    for (i=0;i<M;i++) {
      unew[(i + 1)*N_LEN+N] = unew[(i + 1)*N_LEN];
      vnew[i*N_LEN] = vnew[i*N_LEN+N];
      pnew[i*N_LEN+N] = pnew[i*N_LEN];
    }

    unew[N] = unew[M*N_LEN];
    vnew[M*N_LEN] = vnew[N];
    pnew[M*N_LEN+N] = pnew[0];

    time = time + dt;

    // Time smoothing and update for next cycle

    if ( ncycle > 1 ) {

      c1 = wtime(); 
      for (i=0;i<M_LEN;i++) {
        for (j=0;j<N_LEN;j++) {
          int idx = i*N_LEN+j;
          uold[idx] = u[idx] + alpha * (unew[idx] - 2. * u[idx] + uold[idx]);
          vold[idx] = v[idx] + alpha * (vnew[idx] - 2. * v[idx] + vold[idx]);
          pold[idx] = p[idx] + alpha * (pnew[idx] - 2. * p[idx] + pold[idx]);
        }
      }

      // Dependency

#ifdef _COPY_
      for (i=0;i<M_LEN;i++) {
        for (j=0;j<N_LEN;j++) {
          int idx = i*N_LEN+j;
          u[idx] = unew[idx];
          v[idx] = vnew[idx];
          p[idx] = pnew[idx];
        }
      }
#else
      dswap( (double **)&u, (double **)&unew);
      dswap( (double **)&v, (double **)&vnew);
      dswap( (double **)&p, (double **)&pnew);

#endif

      c2 = wtime(); 
      t300 = t300 + (c2 - c1);
     
    } else {
      tdt = tdt + tdt;

      for (i=0;i<M_LEN;i++) {
        for (j=0;j<N_LEN;j++) {
          int idx = i*N_LEN+j;
          uold[idx] = u[idx];
          vold[idx] = v[idx];
          pold[idx] = p[idx];
        }
      }
      dswap(&unew, &u);
      dswap(&vnew, &v);
      dswap(&pnew, &p);
    }
  } // ** End of time loop ** 

  // have to swap values back for printing

  dswap( (double **)&u, (double **)&unew);
  dswap( (double **)&v, (double **)&vnew);
  dswap( (double **)&p, (double **)&pnew);

  // Output p, u, v fields and run times.
  if(VAL_OUT){
    write_to_file(pnew, M_LEN, N_LEN, "p.bin");
    write_to_file(unew, M_LEN, N_LEN, "u.bin");
    write_to_file(vnew, M_LEN, N_LEN, "v.bin");
    print_to_file(pnew, M_LEN, N_LEN, "p.txt");
    print_to_file(unew, M_LEN, N_LEN, "u.txt");
    print_to_file(vnew, M_LEN, N_LEN, "v.txt");
  }

  if (L_OUT) {
    ptime = time / 3600.;

    printf(" cycle number %d model time in hours %f\n", ITMAX, ptime);
    printf(" diagonal elements of p\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",pnew[i*N_LEN+i]);
    }
    printf("\n diagonal elements of u\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",unew[i*N_LEN+i]);
    }
    printf("\n diagonal elements of v\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",vnew[i*N_LEN+i]);
    }
    printf("\n");

  

    mfs100 = 0.0;
    mfs200 = 0.0;
    mfs300 = 0.0;
    // gdr t100 etc. now an accumulation of all l100 time
    if ( t100 > 0 ) { mfs100 = ITMAX * 24. * M * N / t100 / 1000000; }
    if ( t200 > 0 ) { mfs200 = ITMAX * 26. * M * N / t200 / 1000000; }
    if ( t300 > 0 ) { mfs300 = ITMAX * 15. * M * N / t300 / 1000000; }

    c2 = wtime(); 
    ctime = c2 - tstart;
    tcyc = ctime / ITMAX;

    printf(" cycle number %d total computer time %f time per cycle %f\n", ITMAX, ctime, tcyc);
    printf(" time and megaflops for loop 100 %.6f %.6f\n", t100, mfs100);
    printf(" time and megaflops for loop 200 %.6f %.6f\n", t200, mfs200);
    printf(" time and megaflops for loop 300 %.6f %.6f\n", t300, mfs300);
  }

  free((void *) u);
  free((void *) v);
  free((void *) p);
  free((void *) unew);
  free((void *) vnew);
  free((void *) pnew);
  free((void *) uold);
  free((void *) vold);
  free((void *) pold);
  free((void *) cu);
  free((void *) cv);
  free((void *) z);
  free((void *) h);
  free((void *) psi);

  return(0);
}

void dswap(double **pA, double **pB)
{
  double *pTemp = *pA;
  *pA = *pB;
  *pB = pTemp;
}

void print_to_file(double *array, int tM, int tN, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        return;
    }
    for (int i = 0; i < tM; i++) {
        for (int j = 0; j < tN; j++) {
            fprintf(file, "%f ", array[(i*(tN))+j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void write_to_file(double *array, int tM, int tN, const char *filename) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        return;
    }
    for (int i = 0; i < tM; i++) {
        for (int j = 0; j < tN; j++) {
          fwrite(&array[(i*(tN))+j], sizeof(double), 1, file);
        }
    }
    fclose(file);
}
