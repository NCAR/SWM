
/* Code converted from shallow_base.f90 using F2C-ACC program. 
 * Manually replaced: 
 * - WRITE statements with printf
 * - MOD operator with % 
 * - system_clock with wtime
 * Fixed several of the array references which had x dimension as 1, 
 * instead of M_LEN. 
 * Fixed values set using d and e notation. 
 * (7 June 2011)
#srun ./SWM_gpu | tee test.out
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
#include <string.h>
#ifdef _OPENACC
#include <openacc.h>
#endif 
#ifdef _OPENMP
#include <omp.h>
#endif
#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))

#define TRUE 1
#define FALSE 0
#define M 3072
#define N 3072
#define M_LEN M + 1
#define N_LEN N + 1
#define ITMAX 4000
#define L_OUT TRUE

typedef double real;
//typedef float real;

extern double wtime(); 
extern void dswap(real **a, real **b);
extern void periodic_cont_state_fused(const int m, const int n, real u[m+2][n+2], real v[m+2][n+2], real p[m+2][n+2]);
extern int output_csv_var( char *filename, int m, int n, double* var);
void advance(const int m, const int n, real u[m+2][n+2], real v[m+2][n+2], real p[m+2][n+2],
	     real uold[m+2][n+2], real vold[m+2][n+2], real pold[m+2][n+2],
             real unew[m+2][n+2], real vnew[m+2][n+2], real pnew[m+2][n+2], 
	     real tdt, real dx, real dy);

void update(const int m, const int n, 
	    real u[m+2][n+2], real v[m+2][n+2], real p[m+2][n+2],
	    real uold[m+2][n+2], real vold[m+2][n+2], real pold[m+2][n+2],
            real unew[m+2][n+2], real vnew[m+2][n+2], real pnew[m+2][n+2], 
	    real tdt, real dx, real dy, real alpha);

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
  real u[3][M+2][N+2],v[3][M+2][N+2],p[3][M+2][N+2];
  real psi[M+2][N+2];

#pragma acc enter data copyin(u[:3][:M+2][:N+2],v[:3][:M+2][:N+2],p[:3][:M+2][:N+2],\
                              psi[:M+2][:N+2])
  real dt,tdt,dx,dy,a,alpha,el,pi;
  real tpi,di,dj,pcf;
  real tdts8,tdtsdx,tdtsdy,fsdx,fsdy;

  int mnmin,ncycle;
  int i,j;

    int m;
    int n;
    char id[32] = "";
    if (argc == 2) {
      m = atoi(argv[1]);
      n = atoi(argv[1]);
    }
    else if (argc == 3) {
      m = atoi(argv[1]);
      n = atoi(argv[2]);
    }
    else if (argc == 4) {
      m = atoi(argv[1]);
      n = atoi(argv[2]);
      strcat(id,argv[3]);
    }
    else {
      m = M;
      n = N;
    }

  int new = 2;
  int mid = 1;
  int old = 0;
 
  // timer variables 

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

  el = n * dx;
  pi = 4. * atanf(1.);
  tpi = pi + pi;
  di = tpi / m;
  dj = tpi / n;
  pcf = pi * pi * a * a / (el * el);

#pragma acc data present(u[:3][:m+2][:n+2],v[:3][:m+2][:n+2],p[:3][:m+2][:n+2],\
                              psi[:m+2][:n+2])
#pragma omp parallel for simd private(i,j)
 // Initial values of the stream function and p 
#pragma acc parallel loop present(psi[:m+2][:n+2])
  for (i=0;i<m+1;i++) {
    for (j=0;j<n+1;j++) {
      psi[i][j] = a * sin((i + .5) * di) * sin((j + .5) * dj);
    }
  }
    
  // Initialize velocities
#pragma omp parallel for simd private(i,j)
 #pragma acc parallel loop present(u[:3][:m+2][:n+2],v[:3][:m+2][:n+2],p[:3][:m+2][:n+2])
  for (i=0;i<m;i++) {
    for (j=0;j<n;j++) {
      u[mid][i + 1][j + 1] = -(psi[i + 1][j + 1] - psi[i + 1][j]) / dy;
      v[mid][i + 1][j + 1] = (psi[i + 1][j + 1] - psi[i][j + 1]) / dx;
      p[mid][i + 1][j + 1] = pcf * (cos(2. * (i) * di) + cos(2. * (j) * dj)) + 50000.;
    }
  }
  // Periodic Continuation
  // (in a distributed memory code this would be MPI halo exchanges)
  periodic_cont_state_fused(m, n, u[mid], v[mid], p[mid]);
 
#pragma omp parallel for simd private(i,j)
#pragma acc update host(u[mid:1][:m+2][:n+2],v[mid:1][:m+2][:n+2],p[mid:1][:m+2][:n+2])
#pragma acc parallel loop present(u[:3][:m+2][:n+2],v[:3][:m+2][:n+2],p[:3][:m+2][:n+2]) private(i,j)
  for (i=0;i<m+2;i++) {
    for (j=0;j<n+2;j++) {
      u[old][i][j] = u[mid][i][j];
      v[old][i][j] = v[mid][i][j];
      p[old][i][j] = p[mid][i][j];
    }
  }
     
    // Get difference of p values from 50000
    real dp[(M+2)*(N+2)];
    for (int i=0; i<m+2; i++){
      for (int j=0; j<n+2; j++)
        dp[i*j]=p[mid][i][j]-50000.;
    }

    // Set name of output csv file based on problem size
    char initfile[128] = "swm_init.";
    char size_char[128];
    char tail[128] = "";
    sprintf(size_char, "%d", m);
    strcat(tail,size_char);
    strcat(tail,".");
    sprintf(size_char, "%d", n);
    strcat(tail,size_char);
    strcat(tail,".");
    strcat(tail,id);
    strcat(tail,".csv");
    strcat(initfile,tail);
 
    int outerr = output_csv_var(initfile, m, n, dp);
        
    if (outerr == 0){
    printf("init file output complete\n");
    }

  // Start timer

  time = 0.;  // simulation time

  double tup   = 0.;
  double tcopy = 0.;

  // take first (half) timestep

  advance(m, n, u[mid], v[mid], p[mid], 
                u[old], v[old], p[old], 
                u[new], v[new], p[new], 
                tdt, dx, dy);

  periodic_cont_state_fused(m, n, u[new], v[new], p[new]);

  time = time + tdt;
#pragma omp parallel for simd private(i,j)
#pragma acc update host(u[new:1][:m+2][:n+2],v[new:1][:m+2][:n+2],p[new:1][:m+2][:n+2])
#pragma acc parallel loop collapse(2)\
 present(u[:3][:m+2][:n+2],v[:3][:m+2][:n+2],p[:3][:m+2][:n+2]) private(i,j)
  for (i=0;i<m+2;i++) {
    for (j=0;j<n+2;j++) {
      u[old][i][j] = u[mid][i][j];
      v[old][i][j] = v[mid][i][j];
      p[old][i][j] = p[mid][i][j];
      u[mid][i][j] = u[new][i][j];
      v[mid][i][j] = v[new][i][j];
      p[mid][i][j] = p[new][i][j];
    }
  }

  // Output p, u, v fields after first step.
  if (L_OUT) {
    printf("\n\n");
    printf(" acc diagonal elements of p (1st step)\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",p[new][i+1][i+1]);
    }
    printf("\n acc diagonal elements of u (1st step)\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",u[new][i][i+1]);
    }
    printf("\n acc diagonal elements of v (1st step)\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",v[new][i+1][i]);
    }
    printf("\n\n");
  }

  // From now on, take full timestep 2*dt step

  tdt = tdt + tdt;

  // ** Start of time loop ** 

  tstart = wtime(); 
  for (ncycle=2; ncycle<=ITMAX; ncycle++) {
    
    // Singe update call per step

    c1 = wtime(); 
    update(m, n, u[mid], v[mid], p[mid], 
	   u[old], v[old], p[old], 
	   u[new], v[new], p[new], 
	   tdt, dx, dy, alpha);

    // Periodic Continuation
    // (in a distributed memory code this would be MPI halo exchanges)

    periodic_cont_state_fused(m, n, u[new], v[new], p[new]);
    c2 = wtime();  
    tup = tup + (c2 - c1); 

    // On process mid copy of p[new] halos to p[mid] halos
    // (this is *not* a halo exchange)

    c1 = wtime();  
    // This dependency has to be respected (as currently written)

#define _SWAP_
#ifndef _SWAP_
#pragma omp parallel for simd private(i,j)
#pragma acc parallel loop present(u[:3][:m+2][:n+2],v[:3][:m+2][:n+2],p[:3][:m+2][:n+2])
    for (int i=0;i<m+2;i++) {
      for (int j=0;j<n+2;j++) {
	u[mid][i][j] = u[new][i][j];
	v[mid][i][j] = v[new][i][j];
	p[mid][i][j] = p[new][i][j];
      }
    }
#else
    int tmp;
    tmp = mid;
    mid = new;
    new = tmp;
#endif

    c2 = wtime();  
    tcopy = tcopy + (c2 - c1);

    // update the simulation time

    time = time + dt;

  }  // ** End of timestep loop ** 
     
  // compute timestep loop execution time

  c2 = wtime(); 
  ctime = c2 - tstart;
  tcyc = ctime / (ITMAX-1); //time per cycle

  // Output p, u, v fields and run times.

  if (L_OUT) {
#ifdef _SWAP_
    int tmp;
    tmp = mid;
    mid = new;
    new = tmp;
#endif
#pragma acc update host(u[:3][:m+2][:n+2],v[:3][:m+2][:n+2],p[:3][:m+2][:n+2]) 
    ptime = time / 3600.;
    int nits = ITMAX-1;
    printf(" acc cycle number %d model time in hours %f\n", nits, ptime);
    printf("\n\n");
    printf(" acc diagonal elements of p (%d steps)\n",nits);
    for (i=0; i<mnmin; i++) {
      printf("%f ",p[new][i+1][i+1]);
    }
    printf("\n acc diagonal elements of u (%d steps)\n",nits);
    for (i=0; i<mnmin; i++) {
      printf("%f ",u[new][i][i+1]);
    }
    printf("\n acc diagonal elements of v (%d step)\n",nits);
    for (i=0; i<mnmin; i++) {
      printf("%f ",v[new][i+1][i]);
    }
    printf("\n\n");

    real mflops_up = 0.0;
    real mbps_copy = 0.0;

    // gdr t100 etc. now an accumulation of all l100 time

    if ( tup > 0 )   { mflops_up   = nits * 65. * m * n / tup / 1000000; }
    if ( tcopy > 0 ) { mbps_copy   = nits * sizeof(real)*3*(m+2)*(n+2) / tcopy / 1000000; }
    printf(" acc cycle number %d total computer time %f time per cycle %f\n", nits, ctime, tcyc);
    printf(" acc time and megaflops for update %.6f %.6f\n", tup, mflops_up);
    printf(" acc time and megabytes/sec for copy %.6f %.6f\n", tcopy, mbps_copy);
  }

    // output to .csv file
    
    for (int i=0; i<m+2; i++){
      for (int j=0; j<n+2; j++)
        dp[i*j]=p[new][i][j]-50000.;
    }

    char endfile[32] = "swm_end.";
    strcat(endfile,tail);

    outerr = output_csv_var(endfile, m, n, dp);
    if (outerr == 0){
        printf("end file output complete\n");
    }

  return(0);
}

int output_csv_var( char *filename, int m, int n, double* var  )
{
    FILE *fp;

    fp = fopen(filename, "w+");
    for(int i=1; i<m+1; i++){
       for(int j=1; j<n+1; j++){
           int ij=i*(n+2)+j;
           if (j==n)
               fprintf(fp, "%.15f\n",var[ij]);
           else
               fprintf(fp, "%.15f,",var[ij]);
        }
    }
    fclose(fp);
    return 0;
}

void dswap(real **pA, real **pB)
{
  real *pTemp = *pA;
  *pA = *pB;
  *pB = pTemp;
}

// Periodic continuation of state variables
// Halo structure (x = computed value; o = fill values)
//
// u:
//
// x x x o
// x x x o
// x x x o
// o o o o
//
// v:
//
// o o o o
// o x x x
// o x x x
// o x x x
//
// p:
//
// o o o o
// x x x o
// x x x o
// x x x o
//

void periodic_cont_state_fused(const int m, const int n, real u[m+2][n+2], real v[m+2][n+2], real p[m+2][n+2]){

    int i,j;
    // N/S edge wrapping
#pragma acc enter data copyin(m,n,u[:m+2][:n+2],v[:m+2][:n+2],p[:m+2][:n+2])
#pragma acc  parallel loop  present(u[:m+2][:n+2],v[:m+2][:n+2],p[:m+2][:n+2]) private(j) async 
    for (j=1; j<n+1; j++) {
      u[0  ][j] = u[m][j];
      u[m+1][j] = u[1][j];
      v[0  ][j] = v[m][j];
      v[m+1][j] = v[1][j];
      p[0  ][j] = p[m][j];
      p[m+1][j] = p[1][j];
   }
#pragma acc  parallel loop  present(u[:m+2][:n+2],v[:m+2][:n+2],p[:m+2][:n+2]) private(i) async
    for (i=1; i<m+1; i++) {
#pragma acc cache(u[:][n],u[:][0])
      u[i][0  ] = u[i][n];
      u[i][n+1] = u[i][1];
      v[i][0  ] = v[i][n];
      v[i][n+1] = v[i][1];
      p[i][0  ] = p[i][n];
      p[i][n+1] = p[i][1];
    if(i==m){
      u[0][0]     = u[m][n];
      u[m+1][0]   = u[1][n];
      u[m+1][n+1] = u[1][1];
      u[0][n+1]   = u[m][1];

      v[0][0]     = v[m][n];
      v[m+1][0]   = v[1][n];
      v[m+1][n+1] = v[1][1];
      v[0][n+1]   = v[m][1];

      p[0][0]     = p[m][n];
      p[m+1][0]   = p[1][n];
      p[m+1][n+1] = p[1][1];
      p[0][n+1]   = p[m][1];
     }
   }
}

  // update():
  // Update the state after the first timestep.

  // Compute capital u, capital v, z and h

  // cu values (original):
  // x x x x 0
  // x x x x 0
  // x x x x 0
  // x x x x 0
  // 0 0 0 0 0

  // cu values (fused):
  // x x x x x 
  // x x x x x 
  // x x x x x
  // x x x x x
  // x x x x x
  
  // cv values (original):
  // 0 0 0 0 0
  // 0 x x x x 
  // 0 x x x x 
  // 0 x x x x 
  // 0 x x x x 

  // cv values (fused):
  // x x x x x 
  // x x x x x 
  // x x x x x
  // x x x x x
  // x x x x x

  // z values (original):
  // 0 x x x x 
  // 0 x x x x 
  // 0 x x x x 
  // 0 x x x x 
  // 0 0 0 0 0

  // z values (fused):
  // x x x x x 
  // x x x x x 
  // x x x x x
  // x x x x x
  // x x x x x

  // h values
  // 0 0 0 0 0
  // x x x x 0
  // x x x x 0
  // x x x x 0
  // x x x x 0

  // z values (fused):
  // x x x x x 
  // x x x x x 
  // x x x x x
  // x x x x x
  // x x x x x

  // u values (fused):
  // o o o o o o
  // o x x x x o
  // o x x x x o
  // o x x x x o
  // o o o o o o

  // v values (fused):
  // o o o o o o
  // o x x x x o
  // o x x x x o
  // o x x x x o
  // o x x x x o
  // o o o o o o
  
  // p values (fused):
  // o o o o o o
  // o x x x x o
  // o x x x x o
  // o x x x x o
  // o x x x x o
  // o o o o o o

void update(const int m, const int n, 
	    real u[m+2][n+2], real v[m+2][n+2], real p[m+2][n+2],
	    real uold[m+2][n+2], real vold[m+2][n+2], real pold[m+2][n+2],
            real unew[m+2][n+2], real vnew[m+2][n+2], real pnew[m+2][n+2], 
	    real tdt, real dx, real dy, real alpha){

  int i,j;
#pragma acc data present(u[:m+2][:n+2],v[:m+2][:n+2],p[:m+2][:n+2],\
                        uold[:m+2][:n+2],vold[:m+2][:n+2],pold[:m+2][:n+2],\
                        unew[:m+2][:n+2],vnew[:m+2][:n+2],pnew[:m+2][:n+2])
{
  real fsdx = 4./dx;
  real fsdy = 4./dy;

  // Compute new values u,v and p from cu, cv, z, h

  real tdts8 = tdt / 8.;
  real tdtsdx = tdt / dx;
  real tdtsdy = tdt / dy;

#pragma omp parallel
{
#pragma omp for simd private(i,j)
#pragma acc parallel loop collapse(2) present(u[:m+2][:n+2],v[:m+2][:n+2],p[:m+2][:n+2],\
                             unew[:m+2][:n+2],vnew[:m+2][:n+2],pnew[:m+2][:n+2],\
                             uold[:m+2][:n+2],vold[:m+2][:n+2],pold[:m+2][:n+2]) private(i,j)
  for (int i=1;i<m+1;i++) {
    for (int j=1;j<n+1;j++) {
        int im1 = i-1;
        int jm1 = j-1;
        int ip1 = i+1;
        int jp1 = j+1;
        real p00   = p[i][j];      //10 mem access
        real pp0p1 = p[i][jp1];    //6
        real pp1p0 = p[ip1][j];    //6
        real pp1p1 = p[ip1][jp1];  //3
        real u00   = u[i][j];      //9
        real up0p1  = u[i][jp1];   //4
        real um1p1 = u[im1][jp1];  //4
        real um1p0 = u[im1][j];    //4
        real v00   = v[i][j];      //9
        real vp1p0 = v[ip1][j];    //4
        real vp1m1 = v[ip1][jm1];  //4
        real vp0m1 = v[i][jm1];    //4

        real cut00 = .5 * (pp1p1 + pp0p1)       * up0p1;
        real cut10 = .5 * (pp0p1 + p[im1][jp1]) * um1p1;
        real cut01 = .5 * (pp1p0 + p00)         * u00;
        real cut11 = .5 * (p00   + p[im1][j])   * um1p0;

        real cvt00 = .5 * (pp1p1 + pp1p0)       * vp1p0;
        real cvt10 = .5 * (pp0p1 + p00)         * v00;
        real cvt01 = .5 * (pp1p0 + p[ip1][jm1]) * vp1m1;
        real cvt11 = .5 * (p00   + p[i][jm1])   * vp0m1;

        real zt00 = (fsdx * (vp1p0 - v00)       - fsdy * (up0p1 - u00))      /(p00       + pp1p0       + pp1p1 + pp0p1);
        real zt10 = (fsdx * (v00   - v[im1][j]) - fsdy * (um1p1 - um1p0))    /(p[im1][j] + p00         + pp0p1 + p[im1][jp1]);
        real zt01 = (fsdx * (vp1m1 - vp0m1)     - fsdy * (u00   - u[i][jm1]))/(p[i][jm1] + p[ip1][jm1] + pp1p0 + p00);

        real ht10 = pp0p1 + .25 * (up0p1     * up0p1     + um1p1 * um1p1 + v[i][jp1] * v[i][jp1]  + v00   * v00);
        real ht01 = pp1p0 + .25 * (u[ip1][j] * u[ip1][j] + u00   * u00   + vp1p0     * vp1p0      + vp1m1 * vp1m1);
        real ht11 = p00   + .25 * (u00       * u00       + um1p0 * um1p0 + v00       * v00        + vp0m1 * vp0m1);

        unew[i][j] = uold[i][j] + tdts8  * (zt00  + zt01) * (cvt00 + cvt10 + cvt11 + cvt01) - tdtsdx * (ht01 - ht11);
        vnew[i][j] = vold[i][j] - tdts8  * (zt00  + zt10) * (cut00 + cut10 + cut11 + cut01) - tdtsdy * (ht10 - ht11);
        pnew[i][j] = pold[i][j] - tdtsdx * (cut01 - cut11) - tdtsdy * (cvt10 - cvt11);

	uold[i][j] = u00 + alpha * (unew[i][j] - 2. * u00 + uold[i][j]);
	vold[i][j] = v00 + alpha * (vnew[i][j] - 2. * v00 + vold[i][j]);
	pold[i][j] = p00 + alpha * (pnew[i][j] - 2. * p00 + pold[i][j]);
    }
  }
 } //end omp parallel
} //end data
//#pragma acc exit data delete(cu[:m+1][:n+1], cv[:m+1][:n+1],z[:m+1][:n+1], h[:m+1][:n+1])
}

  // advance: Take the first timestep 
  // Compute capital u, capital v, z and h

  // cu values (original):
  // x x x x 0
  // x x x x 0
  // x x x x 0
  // x x x x 0
  // 0 0 0 0 0

  // cu values (fused):
  // x x x x x 
  // x x x x x 
  // x x x x x
  // x x x x x
  // x x x x x
  
  // cv values (original):
  // 0 0 0 0 0
  // 0 x x x x 
  // 0 x x x x 
  // 0 x x x x 
  // 0 x x x x 

  // cv values (fused):
  // x x x x x 
  // x x x x x 
  // x x x x x
  // x x x x x
  // x x x x x

  // z values (original):
  // 0 x x x x 
  // 0 x x x x 
  // 0 x x x x 
  // 0 x x x x 
  // 0 0 0 0 0

  // z values (fused):
  // x x x x x 
  // x x x x x 
  // x x x x x
  // x x x x x
  // x x x x x

  // h values
  // 0 0 0 0 0
  // x x x x 0
  // x x x x 0
  // x x x x 0
  // x x x x 0

  // z values (fused):
  // x x x x x 
  // x x x x x 
  // x x x x x
  // x x x x x
  // x x x x x

  // u values (fused):
  // o o o o o o
  // o x x x x o
  // o x x x x o
  // o x x x x o
  // o o o o o o

  // v values (fused):
  // o o o o o o
  // o x x x x o
  // o x x x x o
  // o x x x x o
  // o x x x x o
  // o o o o o o
  
  // p values (fused):
  // o o o o o o
  // o x x x x o
  // o x x x x o
  // o x x x x o
  // o x x x x o
  // o o o o o o

void advance(const int m, const int n, 
	     real u[m+2][n+2], real v[m+2][n+2], real p[m+2][n+2],
	     real uold[m+2][n+2], real vold[m+2][n+2], real pold[m+2][n+2],
             real unew[m+2][n+2], real vnew[m+2][n+2], real pnew[m+2][n+2], 
	     real tdt, real dx, real dy){

  real fsdx = 4./dx;
  real fsdy = 4./dy;

  // Compute new values u,v and p from cu, cv, z, h

  real tdts8 = tdt / 8.;
  real tdtsdx = tdt / dx;
  real tdtsdy = tdt / dy;
  int i,j;
#pragma acc enter data copyin(u[:m+2][:n+2],v[:m+2][:n+2],p[:m+2][:n+2],\
                        uold[:m+2][:n+2],vold[:m+2][:n+2],pold[:m+2][:n+2],\
                        unew[:m+2][:n+2],vnew[:m+2][:n+2],pnew[:m+2][:n+2])
#pragma omp parallel
{
#pragma omp for simd private(i,j)
#pragma acc parallel loop collapse(2) present(u[:m+2][:n+2],v[:m+2][:n+2],p[:m+2][:n+2],\
                             unew[:m+2][:n+2],vnew[:m+2][:n+2],pnew[:m+2][:n+2],\
                             uold[:m+2][:n+2],vold[:m+2][:n+2],pold[:m+2][:n+2]) private(i,j)
  for (i=1;i<m+1;i++) {
    for (j=1;j<n+1;j++) {
        int im1 = i-1;
        int jm1 = j-1;
        int ip1 = i+1;
        int jp1 = j+1;
        real p00   = p[i][j];      //10 mem access
        real pp0p1 = p[i][jp1];    //6
        real pp1p0 = p[ip1][j];    //6
        real pp1p1 = p[ip1][jp1];  //3
        real u00   = u[i][j];      //9
        real up0p1  = u[i][jp1];   //4
        real um1p1 = u[im1][jp1];  //4
        real um1p0 = u[im1][j];    //4
        real v00   = v[i][j];      //9
        real vp1p0 = v[ip1][j];    //4
        real vp1m1 = v[ip1][jm1];  //4
        real vp0m1 = v[i][jm1];    //4
        
        real cut00 = .5 * (pp1p1 + pp0p1)       * up0p1;
        real cut10 = .5 * (pp0p1 + p[im1][jp1]) * um1p1;
        real cut01 = .5 * (pp1p0 + p00)         * u00; 
        real cut11 = .5 * (p00   + p[im1][j])   * um1p0;
        
        real cvt00 = .5 * (pp1p1 + pp1p0)       * vp1p0;
        real cvt10 = .5 * (pp0p1 + p00)         * v00;
        real cvt01 = .5 * (pp1p0 + p[ip1][jm1]) * vp1m1;
        real cvt11 = .5 * (p00   + p[i][jm1])   * vp0m1;
        
        real zt00 = (fsdx * (vp1p0 - v00)       - fsdy * (up0p1 - u00))      /(p00       + pp1p0       + pp1p1 + pp0p1);
        real zt10 = (fsdx * (v00   - v[im1][j]) - fsdy * (um1p1 - um1p0))    /(p[im1][j] + p00         + pp0p1 + p[im1][jp1]);
        real zt01 = (fsdx * (vp1m1 - vp0m1)     - fsdy * (u00   - u[i][jm1]))/(p[i][jm1] + p[ip1][jm1] + pp1p0 + p00);
        
        real ht10 = pp0p1 + .25 * (up0p1     * up0p1     + um1p1 * um1p1 + v[i][jp1] * v[i][jp1]  + v00   * v00);
        real ht01 = pp1p0 + .25 * (u[ip1][j] * u[ip1][j] + u00   * u00   + vp1p0     * vp1p0      + vp1m1 * vp1m1);
        real ht11 = p00   + .25 * (u00       * u00       + um1p0 * um1p0 + v00       * v00        + vp0m1 * vp0m1);
        
        unew[i][j] = uold[i][j] + tdts8  * (zt00  + zt01) * (cvt00 + cvt10 + cvt11 + cvt01) - tdtsdx * (ht01 - ht11);
        vnew[i][j] = vold[i][j] - tdts8  * (zt00  + zt10) * (cut00 + cut10 + cut11 + cut01) - tdtsdy * (ht10 - ht11);
        pnew[i][j] = pold[i][j] - tdtsdx * (cut01 - cut11) - tdtsdy * (cvt10 - cvt11); 
        
    }
  }
 } //end omp parallel
}




