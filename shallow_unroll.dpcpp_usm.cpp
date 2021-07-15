// Code converted from shallow_base.f90 using F2C-ACC program
// Manually replaced:
// - WRITE statements with printf
// - MOD operator with %
// - system_clock with wtime
// Fixed several of the array references which had x dimension as 1,
// instead of M_LEN.
// Fixed values set using d and e notation.
// (7 June 2011)

//****************
// 'Pure' C version developed by G.D Riley (UoM) (25 Jan 2012)
// removed all ftocmacros
// used sin and cos not sinf and cosf (since all data are doubles)
// needed to declare arrays +1 to cope with Fortran indexing
//
// Compute on the fly version developed by B. James and C. Miller 18 Dec 2020
//
// Converted to C++ & SYCL by R.D. Loft (5 Feb 2021)
// Starting from hough_transform.cpp Intel OneAPI lab3 code/*
//==============================================================
// hough_transform.cpp:
// Copyright © 2020 Intel Corporation
//
// SPDX-License-Identifier: MIT
// =============================================================

#include <chrono>   // todo: timer --> try high resolution clock to get better results
#include <fstream>  // file
#include <iostream> // io
#include <cmath>
#include <CL/sycl.hpp>
#if FPGA || FPGA_EMULATOR
#include <CL/sycl/INTEL/fpga_extensions.hpp>
#endif

#define TRUE 1
#define FALSE 0
#define M 64
#define N 128
#define M_LEN (M + 2)
#define N_LEN (N + 2)
constexpr size_t  DOMAIN_SIZE = M_LEN*N_LEN;
#define ITMAX 4000
#define L_OUT TRUE

using namespace sycl;


#if 0
// This function reads in a bitmap and outputs an array of pixels
void read_image(char *image_array);

class Hough_Transform_kernel;
#endif

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
//! - Use 8-byte doubles.
//!
//! Further changes (Annette Osprey & Graham Riley, February 2011):
//! - Remove unnecessary halo updates.
//! - Split loops to improve TLB access.
//! - Modify timers to accumulate loop times over whole run.
//! - Remove old-style indentation.
//!
//! Minimal serial version (26 May 2011)


extern void periodic_cont(queue q, int m, int n, double u[DOMAIN_SIZE], double v[DOMAIN_SIZE], double p[DOMAIN_SIZE]);

extern int output_csv_var( char *filename, int m, int n, double* var);
extern void first_step(queue q, int m, int n,
                       double umid[DOMAIN_SIZE], double vmid[DOMAIN_SIZE], double pmid[DOMAIN_SIZE],
                       double uold[DOMAIN_SIZE], double vold[DOMAIN_SIZE], double pold[DOMAIN_SIZE],
                       double unew[DOMAIN_SIZE], double vnew[DOMAIN_SIZE], double pnew[DOMAIN_SIZE],
                       double fsdx, double fsdy, double tdts8, double tdtsdx, double tdtsdy);
void update(queue q, const int m, const int n,
            double u[DOMAIN_SIZE], double v[DOMAIN_SIZE], double p[DOMAIN_SIZE],
            double uold[DOMAIN_SIZE], double vold[DOMAIN_SIZE], double pold[DOMAIN_SIZE],
            double unew[DOMAIN_SIZE], double vnew[DOMAIN_SIZE], double pnew[DOMAIN_SIZE],
            double fsdx, double fsdy, double tdts8, double tdtsdx, double tdtsdy, double alpha);

extern void adv_nsteps(int m, int n,
                       int tlold, int tlmid, int tlnew,
                       double u[3][DOMAIN_SIZE], double v[3][DOMAIN_SIZE], double p[3][DOMAIN_SIZE],
                       double fsdx, double fsdy, double tdts8, double tdtsdx, double tdtsdy,
                       double alpha, int ncycles, double tdt, double time);

void update_time(queue q,
                 double u_in[DOMAIN_SIZE], double v_in[DOMAIN_SIZE], double p_in[DOMAIN_SIZE],
                 double u_out[DOMAIN_SIZE], double v_out[DOMAIN_SIZE], double p_out[DOMAIN_SIZE]);

void init_stream(queue q, int m, int n, double psi[DOMAIN_SIZE]);

void init_velocity(queue q, const int m, const int n, double psi[DOMAIN_SIZE],
                   double u[DOMAIN_SIZE], double v[DOMAIN_SIZE], double p[DOMAIN_SIZE]);

extern void CopyVector(queue &q, double src[DOMAIN_SIZE], double dest[DOMAIN_SIZE]);

extern double wtime();

int main() {
    // Define range, device, and queue
    auto R = range<1>{DOMAIN_SIZE};
    cpu_selector c_selector;
    default_selector d_selector;
    queue q_c(c_selector);
    queue q(d_selector); 
    std::cout << "Device: " << q.get_device().get_info<info::device::name>() << std::endl;

    //Declare state arrays (3 time levels, (M+2)x(N+2) points
    // -- Allocate memory on the host
    double *psi = malloc_host<double>(DOMAIN_SIZE, q_c);

    double **u = malloc_host<double *>(3*DOMAIN_SIZE, q_c);
    double **v = malloc_host<double *>(3*DOMAIN_SIZE, q_c);  
    double **p = malloc_host<double *>(3*DOMAIN_SIZE, q_c);
    for (int i=0; i<3; i++) {
        u[i] = malloc_host<double>(DOMAIN_SIZE, q_c);
        v[i] = malloc_host<double>(DOMAIN_SIZE, q_c);
        p[i] = malloc_host<double>(DOMAIN_SIZE, q_c);
    }
    
    // ** Initialisations **

    int m = M;
    int n = N;

    // Time level indices
    
    int tlnew = 2;
    int tlmid = 1;
    int tlold = 0;
 
    // Note below that two delta t (tdt) is set to dt on the first
    // cycle after which it is reset to dt+dt.
    
    double dt = 90.;
    double tdt = dt;
   
    double dx = 100000.;
    double dy = 100000.;
    double fsdx = 4. / dx;
    double fsdy = 4. / dy;

    double a = 1000000.;
    double alpha = .001;

    double el = n * dx;
    double pi = 4. * atanf(1.);
    double tpi = pi + pi;
    double di = tpi / m;
    double dj = tpi / n;
    double pcf = pi * pi * a * a / (el * el);
    
    double time; // Model time

    std::cout << "stop 0" << std::endl;
    
 // Initial values of the stream function and p
 init_stream(q_c, m, n, psi);
    
 // Initialize velocities
 init_velocity(q_c, m, n, psi, u[tlmid], v[tlmid], p[tlmid]);
     
// Periodic Continuation
// (in a distributed memory code this would be MPI halo exchanges)
periodic_cont(q, m, n, u[tlmid], v[tlmid], p[tlmid]);
std::cout << "stop 1" << std::endl; 
    
update_time(q, u[tlmid], v[tlmid], p[tlmid], u[tlold], v[tlold], p[tlold]);
std::cout << "stop 2" << std::endl;
    
double* dp = new double [DOMAIN_SIZE];
for (int i=0; i<DOMAIN_SIZE; i++){
    dp[i]=p[tlmid][i]-50000.;
    }
char initfile[32] = "swm_init.csv";
 
int outerr = output_csv_var(initfile, m, n, dp);
    
if (outerr == 0){
   std::cout << "init file output complete" << std::endl;
}

double tdts8 = tdt / 8.;
double tdtsdx = tdt / dx;
double tdtsdy = tdt / dy;

double b4first_step = wtime();
  std::cout << "stop 3" << std::endl;  
first_step(q, m, n,
           u[tlmid], v[tlmid], p[tlmid],
           u[tlold], v[tlold], p[tlold],
           u[tlnew], v[tlnew], p[tlnew],
           fsdx, fsdy, tdts8, tdtsdx, tdtsdy);
    std::cout << "stop 4" << std::endl;
 double after_first_step = wtime(); 
std::cout << "first step time: " << after_first_step-b4first_step << std::endl;
    
periodic_cont(q, m, n, u[tlnew], v[tlnew], p[tlnew]);
  
double after_halo = wtime();  
std::cout << "halo update time: " << after_halo-after_first_step << std::endl;
    
time = time + tdt;
update_time(q, u[tlmid], v[tlmid], p[tlmid], u[tlold], v[tlold], p[tlold]);
update_time(q, u[tlnew], v[tlnew], p[tlnew], u[tlmid], v[tlmid], p[tlmid]);
// From now on, take full timestep 2*dt step
tdt = tdt + tdt;
    
tdts8 = tdt / 8.;
tdtsdx = tdt / dx;
tdtsdy = tdt / dy;

// ** Start of time loop **
    
double tup = 0.0;
double tpc = 0.0;
double tstart = wtime();
    
    std::cout << "stop 5" << std::endl;
    
for (int ncycle=2; ncycle<=ITMAX; ncycle++){
        
    // Take a time step
    double c1 = wtime();
    update(q, m, n,
           u[tlmid], v[tlmid], p[tlmid],
           u[tlold], v[tlold], p[tlold],
           u[tlnew], v[tlnew], p[tlnew],
           fsdx, fsdy, tdts8, tdtsdx, tdtsdy, alpha);
    double c2 = wtime();
    //std::cout << "update time: " << c2-c1 << std::endl;
    tup = tup + (c2 - c1);
    
    // Perform periodic continuation
    
    c1 = wtime();
    periodic_cont(q, m, n, u[tlnew], v[tlnew], p[tlnew]);
    c2 = wtime();
    tpc = tpc + (c2-c1);
    
    // update the time level pointers
            
    int tmp = tlmid;
    tlmid = tlnew;
    tlnew = tmp;

    // update the simulation time

    time = time + dt;
        
    }
    std::cout << "stop 6" << std::endl;
double tend = wtime();
double total_time = tend - tstart;
double tcyc = total_time / (ITMAX-1); //time per cycle
    
std::cout << "run complete..." << std::endl;
std::cout << "model time (sec) = " << time << std::endl;

int tmp;
tmp = tlmid;
tlmid = tlnew;
tlnew = tmp;
    
if (L_OUT) {
// if on the GPU, synch data with host
    double ptime = time / 3600.;
    int nits = ITMAX-1;
    int mnmin = std::min(m,n);
    printf(" gpu cycle number %d model time in hours %f\n", nits, ptime);
    printf("\n\n");
    printf(" gpu diagonal elements of p (%d steps)\n",nits+1);
    for (int i=0; i<mnmin; i++) {
      printf("%f ",p[tlnew][(i+1)*(n+2)+(i+1)]);
    }
    printf("\n gpu diagonal elements of u (%d steps)\n",nits+1);
    for (int i=0; i<mnmin; i++) {
      printf("%f ",u[tlnew][i*(n+2)+(i+1)]);
    }
    printf("\n gpu diagonal elements of v (%d step)\n",nits+1);
    for (int i=0; i<mnmin; i++) {
      printf("%f ",v[tlnew][(i+1)*(n+2)+i]);
    }
    printf("\n\n");

    double mflops_up;
    double mbps_pc;

    // gdr t100 etc. now an accumulation of all l100 time

    if ( tup > 0 )   { mflops_up   = nits * 65. * m * n / tup / 1000000; }
    if ( tpc > 0 ) { mbps_pc   = nits * sizeof(double)*3*(2*m+2*n+4) / tpc / 1000000; }
    printf(" cycle number %d total computer time %f, time per cycle %f\n", nits, total_time, tcyc);
    printf(" time (usec/step) %.6f and mflops %.6f\n", 1000000.*tup/nits, mflops_up);
    printf(" time (usec/step) %.6f and mbps %.6f\n", 1000000.*tpc/nits, mbps_pc);
}

// output to .csv file
  
for (int i=0; i<DOMAIN_SIZE; i++){
    dp[i]=p[tlnew][i]-50000.;
    }

char endfile[32] = "swm_h100_dpcpp.csv";
outerr = output_csv_var(endfile, m, n, dp);
if (outerr == 0){
   std::cout << "end file output complete" << std::endl;
   }
    
    // Cleanup
    free(u, q_c);
    free(v, q_c);
    free(p, q_c);
    free(psi, q_c);
    
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
               fprintf(fp, "%.15f;",var[ij]);
           else
               fprintf(fp, "%.15f,",var[ij]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    return 0;
}

void periodic_cont(queue q, const int m, const int n, double u[DOMAIN_SIZE], 
                   double v[DOMAIN_SIZE], double p[DOMAIN_SIZE]) {

    auto R = range<1>{DOMAIN_SIZE};
    
    // -- Allocate memory on the device
    double *u_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *v_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *p_d = malloc_device<double>(DOMAIN_SIZE, q);
    
    // Copy host to device
    CopyVector(q, u, u_d);
    CopyVector(q, v, v_d);
    CopyVector(q, p, p_d);
    
    q.wait(); // Sync
  
        // Corner point exchange indices
        int hsw = 0;
        int hne = (m+1)*(n+2)+n+1;
        int hse = n+1;
        int hnw = (m+1)*(n+2);
    
        int ine = m*(n+2)+n;
        int ise = 1*(n+2)+n;
        int isw = 1*(n+2)+1;
        int inw = m*(n+2)+1;

        q.parallel_for(R, [=](auto index) {
            int j = index%(N+2);
            int i = (int) (index - j)/(N+2);
            if (i==0 || j==0 || i == M+1 || j== N+1) {}
            else {
                // North-South periodic continuation
                // Labeling conventions:
                // h = halo; i = interior;
                // north = m; south = 1;
                int hnorth = (m+1)*(n+2) + j;
                int hsouth = j;
                int inorth = m*(n+2) + j;
                int isouth = (n+2) + j;
        
                u_d[hsouth] = u_d[inorth];
                u_d[hnorth] = u_d[isouth];
        
                v_d[hsouth] = v_d[inorth];
                v_d[hnorth] = v_d[isouth];
        
                p_d[hsouth] = p_d[inorth];
                p_d[hnorth] = p_d[isouth];

                // East-West periodic continuation
                // h = halo; i = interior;
                // east = n; west = 1;
                int hwest = i*(n+2);
                int heast = i*(n+2) + n+1;
                int iwest = i*(n+2) + 1;
                int ieast = i*(n+2) + n;
        
                u_d[hwest] = u_d[ieast];
                u_d[heast] = u_d[iwest];
        
                v_d[hwest] = v_d[ieast];
                v_d[heast] = v_d[iwest];
        
                p_d[hwest] = p_d[ieast];
                p_d[heast] = p_d[iwest];
        
                if (i == m) {
                    u_d[hsw]   = u_d[ine];
                    u_d[hnw]   = u_d[ise];
                    u_d[hne]   = u_d[isw];
                    u_d[hse]   = u_d[inw];

                    v_d[hsw]   = v_d[ine];
                    v_d[hnw]   = v_d[ise];
                    v_d[hne]   = v_d[isw];
                    v_d[hse]   = v_d[inw];
            
                    p_d[hsw]   = p_d[ine];
                    p_d[hnw]   = p_d[ise];
                    p_d[hne]   = p_d[isw];
                    p_d[hse]   = p_d[inw];
                }
            }
        });
    q.wait(); // Synch

    // Copy device to host
    CopyVector(q, u_d, u);
    CopyVector(q, v_d, v);
    CopyVector(q, p_d, p);
    
    q.wait(); // Synch
    
    // Cleanup
    free(u_d, q);
    free(v_d, q);
    free(p_d, q);
}


// numerically the first step is special
void first_step(queue q, const int m, const int n,
         double u[DOMAIN_SIZE], double v[DOMAIN_SIZE], double p[DOMAIN_SIZE],
         double uold[DOMAIN_SIZE], double vold[DOMAIN_SIZE], double pold[DOMAIN_SIZE],
         double unew[DOMAIN_SIZE], double vnew[DOMAIN_SIZE], double pnew[DOMAIN_SIZE],
         double fsdx, double fsdy, double tdts8, double tdtsdx, double tdtsdy) {

    // -- Allocate memory on the device
    double *u_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *v_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *p_d = malloc_device<double>(DOMAIN_SIZE, q);
    
    double *uold_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *vold_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *pold_d = malloc_device<double>(DOMAIN_SIZE, q);
    
    double *unew_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *vnew_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *pnew_d = malloc_device<double>(DOMAIN_SIZE, q);
    
    // Copy host to device
    CopyVector(q, u, u_d);
    CopyVector(q, v, v_d);
    CopyVector(q, p, p_d);
    
    CopyVector(q, uold, uold_d);
    CopyVector(q, vold, vold_d);
    CopyVector(q, pold, pold_d);
    
    CopyVector(q, unew, unew_d);
    CopyVector(q, vnew, vnew_d);
    CopyVector(q, pnew, pnew_d);
    
    q.wait(); // Sync
    
    
    auto R = range<1>{DOMAIN_SIZE};
    q.parallel_for(R, [=](auto ij) {
            int j = ij%(n+2);
            int i = (int) (ij - j)/(n+2);

            // Skip updating the ghost nodes
            if (i==0 || j==0 || i == m+1 || j== n+1) {}
            else {
                int ijp1 = ij+1;
                int ijm1 = ij-1;
                int ip1j = ij+(n+2);
                int ip1jp1 = ip1j+1;
                int ip1jm1 = ip1j-1;
                
                int im1j= ij-(n+2);
                int im1jp1= im1j+1;
                
                double p00   = p_d[ij];      //10 mem access
                double pp0p1 = p_d[ijp1];    //6
                double pp1p0 = p_d[ip1j];    //6
                double pp1p1 = p_d[ip1jp1];  //3
                double u00   = u_d[ij];      //9
                double up0p1 = u_d[ijp1];   //4
                double um1p1 = u_d[im1jp1];  //4
                double um1p0 = u_d[im1j];    //4
                double v00   = v_d[ij];      //9
                double vp1p0 = v_d[ip1j];    //4
                double vp1m1 = v_d[ip1jm1];  //4
                double vp0m1 = v_d[ijm1];    //4
        
                double cut00 = .5 * (pp1p1 + pp0p1)       * up0p1;
                double cut10 = .5 * (pp0p1 + p_d[im1jp1]) * um1p1;
                double cut01 = .5 * (pp1p0 + p00)         * u00;
                double cut11 = .5 * (p00   + p_d[im1j])   * um1p0;
        
                double cvt00 = .5 * (pp1p1 + pp1p0)       * vp1p0;
                double cvt10 = .5 * (pp0p1 + p00)         * v00;
                double cvt01 = .5 * (pp1p0 + p_d[ip1jm1]) * vp1m1;
                double cvt11 = .5 * (p00   + p_d[ijm1])   * vp0m1;
        
                double zt00 = (fsdx * (vp1p0 - v00)     - fsdy * (up0p1 - u00))    /(p00     + pp1p0       + pp1p1 + pp0p1);
                double zt10 = (fsdx * (v00   - v_d[im1j]) - fsdy * (um1p1 - um1p0))  /(p_d[im1j] + p00         + pp0p1 + p_d[im1jp1]);
                double zt01 = (fsdx * (vp1m1 - vp0m1)   - fsdy * (u00   - u_d[ijm1]))/(p_d[ijm1] + p_d[ip1jm1]   + pp1p0 + p00);
        
                double ht10 = pp0p1 + .25 * (up0p1     * up0p1     + um1p1 * um1p1 + v_d[ijp1] * v_d[ijp1]  + v00   * v00);
                double ht01 = pp1p0 + .25 * (u_d[ip1j] * u_d[ip1j] + u00   * u00   + vp1p0     * vp1p0      + vp1m1 * vp1m1);
                double ht11 = p00   + .25 * (u00       * u00       + um1p0 * um1p0 + v00       * v00        + vp0m1 * vp0m1);
        
                unew_d[ij] = uold_d[ij] + tdts8  * (zt00  + zt01) * (cvt00 + cvt10 + cvt11 + cvt01) - tdtsdx * (ht01 - ht11);
                vnew_d[ij] = vold_d[ij] - tdts8  * (zt00  + zt10) * (cut00 + cut10 + cut11 + cut01) - tdtsdy * (ht10 - ht11);
                pnew_d[ij] = pold_d[ij] - tdtsdx * (cut01 - cut11) - tdtsdy * (cvt10 - cvt11);
        
            }
        });
    q.wait();
    
    // Copy device to host
    CopyVector(q, u_d, u);
    CopyVector(q, v_d, v);
    CopyVector(q, p_d, p);
    
    CopyVector(q, uold_d, uold);
    CopyVector(q, vold_d, vold);
    CopyVector(q, pold_d, pold);
    
    CopyVector(q, unew_d, unew);
    CopyVector(q, vnew_d, vnew);
    CopyVector(q, pnew_d, pnew);
    
    q.wait(); // Sync
    
    // Cleanup
    free(u_d, q);
    free(v_d, q);
    free(p_d, q);
    
    free(uold_d, q);
    free(vold_d, q);
    free(pold_d, q);
    
    free(unew_d, q);
    free(vnew_d, q);
    free(pnew_d, q);
}

// update (assumes first step has been called).
void update(queue q, const int m, const int n,
         double u[DOMAIN_SIZE], double v[DOMAIN_SIZE], double p[DOMAIN_SIZE],
         double uold[DOMAIN_SIZE], double vold[DOMAIN_SIZE], double pold[DOMAIN_SIZE],
         double unew[DOMAIN_SIZE], double vnew[DOMAIN_SIZE], double pnew[DOMAIN_SIZE],
         double fsdx, double fsdy, double tdts8, double tdtsdx, double tdtsdy, double alpha) {
    
    // -- Allocate memory on the device
    double *u_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *v_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *p_d = malloc_device<double>(DOMAIN_SIZE, q);
    
    double *uold_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *vold_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *pold_d = malloc_device<double>(DOMAIN_SIZE, q);
    
    double *unew_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *vnew_d = malloc_device<double>(DOMAIN_SIZE, q);
    double *pnew_d = malloc_device<double>(DOMAIN_SIZE, q);
    
    // Copy host to device
    CopyVector(q, u, u_d);
    CopyVector(q, v, v_d);
    CopyVector(q, p, p_d);
    
    CopyVector(q, uold, uold_d);
    CopyVector(q, vold, vold_d);
    CopyVector(q, pold, pold_d);
    
    CopyVector(q, unew, unew_d);
    CopyVector(q, vnew, vnew_d);
    CopyVector(q, pnew, pnew_d);
    
    q.wait(); // Sync
    
    auto R = range<1>{DOMAIN_SIZE};
    q.parallel_for(R, [=](auto ij) {
            int j = ij%(n+2);
            int i = (int) (ij - j)/(n+2);

            // Skip updating the ghost nodes
            if (i==0 || j==0 || i == m+1 || j== n+1) {}
            else {
                int ijp1 = ij+1;
                int ijm1 = ij-1;
        
                int ip1j   = ij+(n+2);
                int ip1jp1 = ip1j+1;
                int ip1jm1 = ip1j-1;
        
                int im1j   = ij-(n+2);
                int im1jp1 = im1j+1;
        
                double p00   = p_d[ij];      //10 mem access
                double pp0p1 = p_d[ijp1];    //6
                double pp1p0 = p_d[ip1j];    //6
                double pp1p1 = p_d[ip1jp1];  //3
                double u00   = u_d[ij];      //9
                double up0p1 = u_d[ijp1];    //4
                double um1p1 = u_d[im1jp1];  //4
                double um1p0 = u_d[im1j];    //4
                double v00   = v_d[ij];      //9
                double vp1p0 = v_d[ip1j];    //4
                double vp1m1 = v_d[ip1jm1];  //4
                double vp0m1 = v_d[ijm1];    //4
        
                double cut00 = .5 * (pp1p1 + pp0p1)     * up0p1;
                double cut10 = .5 * (pp0p1 + p_d[im1jp1]) * um1p1;
                double cut01 = .5 * (pp1p0 + p00)       * u00;
                double cut11 = .5 * (p00   + p_d[im1j])   * um1p0;
        
                double cvt00 = .5 * (pp1p1 + pp1p0)     * vp1p0;
                double cvt10 = .5 * (pp0p1 + p00)       * v00;
                double cvt01 = .5 * (pp1p0 + p_d[ip1jm1]) * vp1m1;
                double cvt11 = .5 * (p00   + p_d[ijm1])   * vp0m1;
        
                double zt00 = (fsdx * (vp1p0 - v00)     - fsdy * (up0p1 - u00))     /(p00     + pp1p0     + pp1p1 + pp0p1);
                double zt10 = (fsdx * (v00   - v_d[im1j]) - fsdy * (um1p1 - um1p0))   /(p_d[im1j] + p00       + pp0p1 + p_d[im1jp1]);
                double zt01 = (fsdx * (vp1m1 - vp0m1)   - fsdy * (u00   - u_d[ijm1])) /(p_d[ijm1] + p_d[ip1jm1] + pp1p0 + p00);
        
                double ht10 = pp0p1 + .25 * (up0p1     * up0p1     + um1p1 * um1p1 + v_d[ijp1] * v_d[ijp1]  + v00   * v00);
                double ht01 = pp1p0 + .25 * (u_d[ip1j]   * u_d[ip1j]   + u00   * u00   + vp1p0   * vp1p0    + vp1m1 * vp1m1);
                double ht11 = p00   + .25 * (u00       * u00       + um1p0 * um1p0 + v00     * v00      + vp0m1 * vp0m1);
        
                unew_d[ij] = uold_d[ij] + tdts8  * (zt00  + zt01) * (cvt00 + cvt10 + cvt11 + cvt01) - tdtsdx * (ht01 - ht11);
                vnew_d[ij] = vold_d[ij] - tdts8  * (zt00  + zt10) * (cut00 + cut10 + cut11 + cut01) - tdtsdy * (ht10 - ht11);
                pnew_d[ij] = pold_d[ij] - tdtsdx * (cut01 - cut11) - tdtsdy * (cvt10 - cvt11);
        
                uold_d[ij] = u00 + alpha * (unew_d[ij]  - 2. * u00 + uold_d[ij]);
                vold_d[ij] = v00 + alpha * (vnew_d[ij]  - 2. * v00 + vold_d[ij]);
                pold_d[ij] = p00 + alpha * (pnew_d[ij]  - 2. * p00 + pold_d[ij]);
            }
        });
    q.wait();
    
    // Copy device to host
    CopyVector(q, u_d, u);
    CopyVector(q, v_d, v);
    CopyVector(q, p_d, p);
    
    CopyVector(q, uold_d, uold);
    CopyVector(q, vold_d, vold);
    CopyVector(q, pold_d, pold);
    
    CopyVector(q, unew_d, unew);
    CopyVector(q, vnew_d, vnew);
    CopyVector(q, pnew_d, pnew);
    
    q.wait(); // Sync
    
    // Cleanup
    free(u_d, q);
    free(v_d, q);
    free(p_d, q);
    
    free(uold_d, q);
    free(vold_d, q);
    free(pold_d, q);
    
    free(unew_d, q);
    free(vnew_d, q);
    free(pnew_d, q);
}

void update_time(queue q, 
                 double u_in[DOMAIN_SIZE], double v_in[DOMAIN_SIZE], double p_in[DOMAIN_SIZE],
                 double u_out[DOMAIN_SIZE], double v_out[DOMAIN_SIZE], double p_out[DOMAIN_SIZE]) {
    
    double *u_d_in = malloc_device<double>(DOMAIN_SIZE, q);
    double *v_d_in = malloc_device<double>(DOMAIN_SIZE, q);
    double *p_d_in = malloc_device<double>(DOMAIN_SIZE, q);
    
    double *u_d_out = malloc_device<double>(DOMAIN_SIZE, q);
    double *v_d_out = malloc_device<double>(DOMAIN_SIZE, q);
    double *p_d_out = malloc_device<double>(DOMAIN_SIZE, q);

    // Copy host to device
    CopyVector(q, u_in, u_d_in);
    CopyVector(q, v_in, v_d_in);
    CopyVector(q, p_in, p_d_in);

    CopyVector(q, u_out, u_d_out);
    CopyVector(q, v_out, v_d_out);
    CopyVector(q, p_out, p_d_out);
    
    q.wait(); // Sync
        
    auto R = range<1>{DOMAIN_SIZE};

    q.parallel_for(R, [=](auto ij) {
            u_d_out[ij] = u_d_in[ij];
            v_d_out[ij] = v_d_in[ij];
            p_d_out[ij] = p_d_in[ij];   
        });
    q.wait(); // Synch

    // Copy device to host
    CopyVector(q, u_d_in, u_in);
    CopyVector(q, v_d_in, v_in);
    CopyVector(q, p_d_in, p_in);

    CopyVector(q, u_d_out, u_out);
    CopyVector(q, v_d_out, v_out);
    CopyVector(q, p_d_out, p_out);
    
    q.wait(); // Synch
    
    // Cleanup
    free(u_d_in, q);
    free(v_d_in, q);
    free(p_d_in, q);
    
    free(u_d_out, q);
    free(v_d_out, q);
    free(p_d_out, q);
}

void init_stream(queue q, int m, int n, double psi[DOMAIN_SIZE]) {
    
    double a = 1000000.;
    double alpha = .001;
    double dx = 100000.;
    double dy = 100000.;
    double el = n * dx;
    double pi = 4. * atanf(1.);
    double tpi = pi + pi;
    double di = tpi / m;
    double dj = tpi / n;
    double pcf = pi * pi * a * a / (el * el);

    auto R = range<1>{DOMAIN_SIZE};

    auto e = q.parallel_for(R, [=](auto ij) {
        int j = ij%(n+2);
        int i = (int) (ij - j)/(n+2);
        
        // Skip updating the ghost nodes
        if (i == m+1 || j== n+1) {}
        else {
            psi[ij] = a * sin((i + .5) * di) * sin((j + .5) * dj);
        }
    });
    e.wait();
}

void init_velocity(queue q, const int m, const int n, double psi[DOMAIN_SIZE],
                   double u[DOMAIN_SIZE], double v[DOMAIN_SIZE], double p[DOMAIN_SIZE]) {
    double dx = 100000.;
    double dy = 100000.;
    double a = 1000000.;
    double alpha = .001;
    double el = n * dx;
    double pi = 4. * atanf(1.);
    double tpi = pi + pi;
    double di = tpi / m;
    double dj = tpi / n;
    double pcf = pi * pi * a * a / (el * el);
    
    auto R = range<1>{DOMAIN_SIZE};
    auto e = q.parallel_for(R, [=](auto ij) {
        int j = ij%(n+2);
        int i = (int) (ij - j)/(n+2);
        
        int ijm1 = ij-1;
        int im1j= ij-(n+2);
        
        if (i==0 || j==0 || i == m+1 || j== n+1) {}
        else {
            u[ij] = -(psi[ij] - psi[ijm1]) / dy;
            v[ij] = (psi[ij] - psi[im1j]) / dx;
            p[ij] = pcf * (cos(2. * (i-1) * di) + cos(2. * (j-1) * dj)) + 50000.;
        }
    });
    e.wait();
}

void CopyVector(queue &q, double src[DOMAIN_SIZE], double dest[DOMAIN_SIZE]) {
    
    q.submit([&](handler &h) {
      h.memcpy(dest, &src[0], DOMAIN_SIZE * sizeof(double)); 
    });
}
