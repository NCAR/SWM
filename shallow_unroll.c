
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
#include <string.h>
#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))

#define TRUE 1
#define FALSE 0
#define M 64
#define N 128
#define M_LEN M + 1
#define N_LEN N + 1
#define ITMAX 4000
#define L_OUT TRUE

typedef double real;
//typedef float real;

extern double wtime(); 
extern void dswap(real **a, real **b);
extern void periodic_cont_state_fused(const int m, const int n, real u[m+2][n+2], real v[m+2][n+2], real p[m+2][n+2]);
extern int output_csv_var( char *filename, int m, int n, real var[m+2][n+2]);
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

    int m;
    int n;
    char id[128] = "";
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

  // solution arrays
  real u[3][m+2][n+2],v[3][m+2][n+2],p[3][m+2][n+2];
  real psi[m+2][n+2];

  real dt,tdt,dx,dy,a,alpha,el,pi;
  real tpi,di,dj,pcf;
  real tdts8,tdtsdx,tdtsdy,fsdx,fsdy;

  int mnmin,ncycle;
  int i,j;

  // Time level indices

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

  // Initial values of the stream function and p
  for (i=0;i<m+1;i++) {
    for (j=0;j<n+1;j++) {
      psi[i][j] = a * sin((i + .5) * di) * sin((j + .5) * dj);
    }
  }
    
  // Initialize velocities
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

  for (i=0;i<m+2;i++) {
    for (j=0;j<n+2;j++) {
      u[old][i][j] = u[mid][i][j];
      v[old][i][j] = v[mid][i][j];
      p[old][i][j] = p[mid][i][j];
    }
  }
     
  // Print initial values
  if ( L_OUT ) {
    printf(" number of points in the x direction %d\n", n); 
    printf(" number of points in the y direction %d\n", m); 
    printf(" grid spacing in the x direction     %f\n", dx); 
    printf(" grid spacing in the y direction     %f\n", dy); 
    printf(" time step                           %f\n", dt); 
    printf(" time filter parameter               %f\n", alpha); 

    mnmin = MIN(m,n);
    printf("\n\n");
    printf(" initial diagonal elements of p\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",p[mid][i+1][i+1]);
    }
    // print a patch of u
#if 0
    printf("\n initial diagonal elements of u\n");
    for (i=0; i<4; i++) {
      for (j=0; j<4; j++) {
	printf("%f ",u[mid][i][j+1]);
      }
      printf("\n");
    }
#endif
    printf("\n initial diagonal elements of u\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",u[mid][i][i+1]);
    }
#if 0
    // print a patch of v
    printf("\n initial diagonal elements of v\n");
    for (i=0; i<4; i++) {
      for (j=0; j<4; j++) {
	printf("%f ",v[mid][i+1][j]);
      }
      printf("\n");
    }
#endif
    printf("\n initial diagonal elements of v\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",v[mid][i+1][i]);
    }
    printf("\n");
  }

    // Get difference of p values from 50000
    real dp[m+2][n+2];
    for (i=0; i<m+2; i++){
      for (j=0; j<n+2; j++)
        dp[i][j]=p[mid][i][j]-50000.;
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
    printf(" diagonal elements of p (1st step)\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",p[new][i+1][i+1]);
    }
    printf("\n diagonal elements of u (1st step)\n");
    for (i=0; i<mnmin; i++) {
      printf("%f ",u[new][i][i+1]);
    }
    printf("\n diagonal elements of v (1st step)\n");
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
    ptime = time / 3600.;
    int nits = ITMAX-1;
    printf(" update steps taken =  %d simulated time in hours %f\n", nits, ptime);
    printf("\n\n");
    printf(" diagonal elements of p (after %d steps)\n",nits+1);
    for (i=0; i<mnmin; i++) {
      printf("%f ",p[new][i+1][i+1]);
    }
    printf("\n diagonal elements of u (after %d steps)\n",nits+1);
    for (i=0; i<mnmin; i++) {
      printf("%f ",u[new][i][i+1]);
    }
    printf("\n diagonal elements of v (after %d step)\n",nits+1);
    for (i=0; i<mnmin; i++) {
      printf("%f ",v[new][i+1][i]);
    }
    printf("\n\n");

    real mflops_up = 0.0;
    real mbps_copy = 0.0;

    // gdr t100 etc. now an accumulation of all l100 time

    if ( tup > 0 )   { mflops_up   = nits * 65. * m * n / tup / 1000000; }
    if ( tcopy > 0 ) { mbps_copy   = nits * sizeof(real)*3*(m+2)*(n+2) / tcopy / 1000000; }
    printf(" cycle number %d total computer time %f time per cycle %f\n", nits, ctime, tcyc);
    printf(" time and megaflops for update %.6f %.6f\n", tup, mflops_up);
    printf(" time and megabytes/sec for copy %.6f %.6f\n", tcopy, mbps_copy);
  }

 // output to .csv file
    
    for (i=0; i<m+2; i++){
      for (j=0; j<n+2; j++)
        dp[i][j]=p[new][i][j]-50000.;
    }

    char endfile[128] = "swm_end.";
    strcat(endfile,tail);

    outerr = output_csv_var(endfile, m, n, dp);
    if (outerr == 0){
        printf("end file output complete\n");
    }

  return(0);
}

int output_csv_var( char *filename, int m, int n, real var[m+2][n+2] )
{
    FILE *fp;

    fp = fopen(filename, "w+");
    int i,j;
    for(i=1; i<m+1; i++){
       for(j=1; j<n+1; j++){
           if (j==n)
               fprintf(fp, "%.15f\n",var[i][j]);
           else
               fprintf(fp, "%.15f,",var[i][j]);
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
    
    for (j=1; j<n+1; j++) {
      u[0  ][j] = u[m][j];
      u[m+1][j] = u[1][j];
      v[0  ][j] = v[m][j];
      v[m+1][j] = v[1][j];
      p[0  ][j] = p[m][j];
      p[m+1][j] = p[1][j];
    }

    // E/W edge wrapping

    for (i=1; i<m+1; i++) {
      u[i][0  ] = u[i][n];
      u[i][n+1] = u[i][1];
      v[i][0  ] = v[i][n];
      v[i][n+1] = v[i][1];
      p[i][0  ] = p[i][n];
      p[i][n+1] = p[i][1];
    }
    
    // Corners (periodic BC's)

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

  real cu[m+1][n+1], cv[m+1][n+1];
  real z[m+1][n+1], h[m+1][n+1];

  real fsdx = 4./dx;
  real fsdy = 4./dy;

  // Compute new values u,v and p from cu, cv, z, h

  real tdts8 = tdt / 8.;
  real tdtsdx = tdt / dx;
  real tdtsdy = tdt / dy;

  int i,j;
  for (i=0;i<m+1;i++) {
    for (j=0;j<n+1;j++) {

      //orig:     cu[i + 1][j] = .5 * (p[i + 1][j] + p[i][j]) * u[i + 1][j];
      cu[i][j] = .5 * (p[i+1][j+1] + p[i][j+1]) * u[i][j+1];

      //orig:      cv[i][j + 1] = .5 * (p[i][j + 1] + p[i][j]) * v[i][j + 1];
      cv[i][j] = .5 * (p[i+1][j+1] + p[i+1][j]) * v[i+1][j];

      //orig:      z[i + 1][j + 1] = (fsdx * (v[i + 1][j + 1] - v[i][j + 1]) - fsdy * (u[i + 1][j + 1] - u[i + 1][j])) / (p[i][j] + p[i + 1][j] + p[i + 1][j + 1] + p[i][j + 1]);
      z[i][j] = (fsdx * (v[i + 1][j] - v[i][j]) - fsdy * (u[i][j + 1] - u[i][j]))/(p[i][j] + p[i + 1][j] + p[i + 1][j + 1] + p[i][j + 1]);

      //orig:      h[i][j] = p[i][j] + .25 * (u[i + 1][j] * u[i + 1][j] + u[i][j] * u[i][j] + v[i][j + 1] * v[i][j + 1] + v[i][j] * v[i][j]);
      h[i][j] = p[i+1][j+1] + .25 * (u[i+1][j+1] * u[i+1][j+1] + u[i][j+1] * u[i][j+1] + v[i+1][j+1] * v[i+1][j+1] + v[i+1][j] * v[i+1][j]);
      //#define _FUSE_LOOPS_
#ifndef _FUSE_LOOPS_
    }
  }

  for (i=1;i<m+1;i++) {
    for (j=1;j<n+1;j++) {
#else
      if (i>0 && j>0){
#endif
	unew[i][j] = uold[i][j] + tdts8 * (z[i][j] + z[i][j-1]) * (cv[i][j] + cv[i-1][j] + cv[i-1][j-1] + cv[i][j-1]) - tdtsdx * (h[i][j-1] - h[i-1][j-1]);
        vnew[i][j] = vold[i][j] - tdts8 * (z[i][j] + z[i-1][j]) * (cu[i][j] + cu[i-1][j] + cu[i-1][j-1] + cu[i][j-1]) - tdtsdy * (h[i-1][j] - h[i-1][j-1]);
	pnew[i][j] = pold[i][j] - tdtsdx * (cu[i][j-1] - cu[i-1][j-1]) - tdtsdy * (cv[i-1][j] - cv[i-1][j-1]); 

	uold[i][j] = u[i][j] + alpha * (unew[i][j] - 2. * u[i][j] + uold[i][j]);
	vold[i][j] = v[i][j] + alpha * (vnew[i][j] - 2. * v[i][j] + vold[i][j]);
	pold[i][j] = p[i][j] + alpha * (pnew[i][j] - 2. * p[i][j] + pold[i][j]);
#ifdef _FUSE_LOOPS
      }
#endif
    }
  }


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

  real cu[m+1][n+1], cv[m+1][n+1];
  real z[m+1][n+1], h[m+1][n+1];

  real fsdx = 4./dx;
  real fsdy = 4./dy;

  // Compute new values u,v and p from cu, cv, z, h

  real tdts8 = tdt / 8.;
  real tdtsdx = tdt / dx;
  real tdtsdy = tdt / dy;

  int i,j;
  for (i=0;i<m+1;i++) {
    for (j=0;j<n+1;j++) {
      //orig:     cu[i + 1][j] = .5 * (p[i + 1][j] + p[i][j]) * u[i + 1][j];
      cu[i][j] = .5 * (p[i+1][j+1] + p[i][j+1]) * u[i][j+1];

      //orig:      cv[i][j + 1] = .5 * (p[i][j + 1] + p[i][j]) * v[i][j + 1];
      cv[i][j] = .5 * (p[i+1][j+1] + p[i+1][j]) * v[i+1][j];

      //orig:      z[i + 1][j + 1] = (fsdx * (v[i + 1][j + 1] - v[i][j + 1]) - fsdy * (u[i + 1][j + 1] - u[i + 1][j])) / (p[i][j] + p[i + 1][j] + p[i + 1][j + 1] + p[i][j + 1]);
      z[i][j] = (fsdx * (v[i + 1][j] - v[i][j]) - fsdy * (u[i][j + 1] - u[i][j]))/(p[i][j] + p[i + 1][j] + p[i + 1][j + 1] + p[i][j + 1]);

      //orig:      h[i][j] = p[i][j] + .25 * (u[i + 1][j] * u[i + 1][j] + u[i][j] * u[i][j] + v[i][j + 1] * v[i][j + 1] + v[i][j] * v[i][j]);
      h[i][j] = p[i+1][j+1] + .25 * (u[i+1][j+1] * u[i+1][j+1] + u[i][j+1] * u[i][j+1] + v[i+1][j+1] * v[i+1][j+1] + v[i+1][j] * v[i+1][j]);

      if (i>0 && j>0){
	unew[i][j] = uold[i][j] + tdts8 * (z[i][j] + z[i][j-1]) * (cv[i][j] + cv[i-1][j] + cv[i-1][j-1] + cv[i][j-1]) - tdtsdx * (h[i][j-1] - h[i-1][j-1]);
        vnew[i][j] = vold[i][j] - tdts8 * (z[i][j] + z[i-1][j]) * (cu[i][j] + cu[i-1][j] + cu[i-1][j-1] + cu[i][j-1]) - tdtsdy * (h[i-1][j] - h[i-1][j-1]);
	pnew[i][j] = pold[i][j] - tdtsdx * (cu[i][j-1] - cu[i-1][j-1]) - tdtsdy * (cv[i-1][j] - cv[i-1][j-1]); 
      }

    }
  }

#ifdef _DEBUG_
  printf("\n advance: cu\n");
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      printf("%.8f ",cu[i][j]);
    }
    printf("\n");
  }

  printf("\n advance: cv\n");
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      printf("%.8f ",cv[i][j]);
    }
    printf("\n");
  }

  printf("\n advance: z\n");
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      printf("%.8e ",z[i][j]);
    }
    printf("\n");
  }

  printf("\n advance: h\n");
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      printf("%.8f ",h[i][j]);
    }
    printf("\n");
  }

  printf("\n advance: unew\n");
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      printf("%.8f ",unew[i+1][j+1]);
    }
    printf("\n");
  }

  printf("\n advance: vnew\n");
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      printf("%.8f ",vnew[i+1][j+1]);
    }
    printf("\n");
  }
  printf("\n advance: uold\n");
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      printf("%.8f ",uold[i+1][j+1]);
    }
    printf("\n");
  }

  printf("\n advance: vold\n");
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      printf("%.8f ",vold[i+1][j+1]);
    }
    printf("\n");
  }
#endif

}




