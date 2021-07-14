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

#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>

#define TRUE 1
#define FALSE 0
const int M = 64;
const int N = 128;
//#define M_LEN (M + 2)
//#define N_LEN (N + 2)
//#define DOMAIN_SIZE M_LEN*N_LEN
#define ITMAX 4000
#define L_OUT TRUE

using namespace std;

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

using namespace std;
#include <cmath>
#include <Kokkos_Core.hpp>

//Create View types which use default execution space
typedef Kokkos::View<double*[3]> AllLevelType;
typedef Kokkos::View<double*> SingleLevelType;

extern void periodic_cont(int m, int n, int tl,
                          AllLevelType uAll, AllLevelType vAll, AllLevelType pAll);

extern int output_csv_var( char *filename, int m, int n, double* var);

extern void first_step(int m, int n, int tlmid, int tlold, int tlnew,
                       AllLevelType uAll, AllLevelType vAll, AllLevelType pAll,
                       double fsdx, double fsdy, double tdts8, double tdtsdx, double tdtsdy);

void update(const int m, const int n, int tlmid, int tlold, int tlnew,
            AllLevelType uAll, AllLevelType vAll, AllLevelType pAll,
            double fsdx, double fsdy, double tdts8, double tdtsdx, double tdtsdy, double alpha);

extern void adv_nsteps(int m, int n, int tlmid, int tlold, int tlnew,
                       AllLevelType uAll, AllLevelType vAll, AllLevelType pAll,
                       double fsdx, double fsdy, double tdts8, double tdtsdx, double tdtsdy,
                       double alpha, int ncycles, double tdt, double time);

extern double wtime();

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    //Print default execution space
    printf("Kokkos execution space: %s\n",
         typeid(Kokkos::DefaultExecutionSpace).name());

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

    cout << "number of points in the x direction " << m << endl;
    cout << "number of points in the y direction " << n << endl;

    const int M_LEN = (m + 2);
    const int N_LEN = (n + 2);
    const int DOMAIN_SIZE = M_LEN * N_LEN;

    Kokkos::Profiling::pushRegion("allocation");
    //Declare state arrays as views ((M+2)x(N+2) points, 3 time levels)
    //Domain size is first dimension to optimize memory layouts
    //Allocate device memory
    AllLevelType u("u",DOMAIN_SIZE);
    AllLevelType v("v",DOMAIN_SIZE);
    AllLevelType p("p",DOMAIN_SIZE);
    SingleLevelType psi("psi",DOMAIN_SIZE);
    // Allocate host memory
    AllLevelType::HostMirror u_h = Kokkos::create_mirror_view(u);
    AllLevelType::HostMirror v_h = Kokkos::create_mirror_view(v);
    AllLevelType::HostMirror p_h = Kokkos::create_mirror_view(p);
    SingleLevelType::HostMirror psi_h = Kokkos::create_mirror_view(psi);
    Kokkos::Profiling::popRegion();
    
    // ** Initialisations **

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
    
    Kokkos::Profiling::pushRegion("initialization");
    // Initial values of the stream function and p
    Kokkos::parallel_for("init_stream_func_vals",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{m+1,n+1}),
      KOKKOS_LAMBDA(const int i, const int j){
        int ii = i*(n+2)+j;
        psi(ii) = a * sin((i + .5) * di) * sin((j + .5) * dj);
      }
    );

    // Initialize velocities
    Kokkos::parallel_for("init_velocities",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{m,n}),
      KOKKOS_LAMBDA(const int i, const int j){
        int ipjp = (i+1)*N_LEN + j+1;
        int ipj = (i+1)*N_LEN + j;
        int ijp = i*N_LEN + j+1;
        u(ipjp,tlmid) = -(psi(ipjp) - psi(ipj)) / dy;
        v(ipjp,tlmid) = (psi(ipjp) - psi(ijp)) / dx;
        p(ipjp,tlmid) = pcf * (cos(2. * (i) * di) + cos(2. * (j) * dj)) + 50000.;
      }
    );
    Kokkos::Profiling::popRegion();
    
    Kokkos::Profiling::pushRegion("setup");
    // Periodic Continuation
    // (in a distributed memory code this would be MPI halo exchanges)
    periodic_cont(m, n, tlmid, u, v, p);
    
    Kokkos::parallel_for("mid_to_old",
      DOMAIN_SIZE,
      KOKKOS_LAMBDA(const int ij){
        u(ij,tlold) = u(ij,tlmid);
        v(ij,tlold) = v(ij,tlmid);
        p(ij,tlold) = p(ij,tlmid);
      }
    );
    Kokkos::deep_copy(p_h,p);
    
    double* dp = new double [DOMAIN_SIZE];
    for (int i=0; i<DOMAIN_SIZE; i++){
      dp[i]=p_h(i,tlmid)-50000.;
    }
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
      std::cout << "init file output complete" << std::endl;
    }

    double tdts8 = tdt / 8.;
    double tdtsdx = tdt / dx;
    double tdtsdy = tdt / dy;
  
    first_step(m, n, tlmid, tlold, tlnew, u, v, p,
               fsdx, fsdy, tdts8, tdtsdx, tdtsdy);
    
    periodic_cont(m, n, tlnew, u, v, p);
    
    time = time + tdt;
    Kokkos::parallel_for("mid-to-old_new-to-mid",
      DOMAIN_SIZE,
      KOKKOS_LAMBDA(const int ij){
        u(ij,tlold) = u(ij,tlmid);
        v(ij,tlold) = v(ij,tlmid);
        p(ij,tlold) = p(ij,tlmid);
        u(ij,tlmid) = u(ij,tlnew);
        v(ij,tlmid) = v(ij,tlnew);
        p(ij,tlmid) = p(ij,tlnew);
      }
    );

    // From now on, take full timestep 2*dt step
    tdt = tdt + tdt;
        
    tdts8 = tdt / 8.;
    tdtsdx = tdt / dx;
    tdtsdy = tdt / dy;

    // ** Start of time loop ** 
    double tup = 0.0;
    double tpc = 0.0;
    double tstart = wtime();
    Kokkos::Profiling::popRegion();
    for (int ncycle=2; ncycle<=ITMAX; ncycle++){
      Kokkos::Profiling::pushRegion("iterations");
      // Take a time step
      double c1 = wtime();
      Kokkos::Profiling::pushRegion("update-function");
      update(m, n, tlmid, tlold, tlnew,
             u, v, p,
             fsdx, fsdy, tdts8, tdtsdx, tdtsdy, alpha);
      Kokkos::Profiling::popRegion();
      double c2 = wtime();
      tup = tup + (c2 - c1);
      
      // Perform periodic continuation
      c1 = wtime();
      Kokkos::Profiling::pushRegion("periodic_cont-function");
      periodic_cont(m, n, tlnew, u, v, p);
      Kokkos::Profiling::popRegion();
      c2 = wtime();
      tpc = tpc + (c2-c1);
      
      // update the time level pointers 
      int tmp = tlmid;
      tlmid = tlnew;
      tlnew = tmp;

      // update the simulation time
      time = time + dt; 
      Kokkos::Profiling::popRegion();
    }
    Kokkos::Profiling::pushRegion("output");
    Kokkos::fence();
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
    // if on the GPU, sync data with host
      Kokkos::deep_copy(p_h,p);
      Kokkos::deep_copy(u_h,u);
      Kokkos::deep_copy(v_h,v);
      double ptime = time / 3600.;
      int nits = ITMAX-1;
      int mnmin = std::min(m,n);
      printf(" gpu cycle number %d model time in hours %f\n", nits, ptime);
      printf("\n\n");
      printf(" gpu diagonal elements of p (%d steps)\n",nits+1);
      for (int i=0; i<mnmin; i++) {
        printf("%f ",p_h((i+1)*(n+2)+(i+1),tlnew));
      }
      printf("\n gpu diagonal elements of u (%d steps)\n",nits+1);
      for (int i=0; i<mnmin; i++) {
        printf("%f ",u_h(i*(n+2)+(i+1),tlnew));
      }
      printf("\n gpu diagonal elements of v (%d step)\n",nits+1);
      for (int i=0; i<mnmin; i++) {
        printf("%f ",v_h((i+1)*(n+2)+i,tlnew));
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
      dp[i]=p_h(i,tlnew)-50000.;
    }

    char endfile[32] = "swm_end.";
    strcat(endfile,tail);

    outerr = output_csv_var(endfile, m, n, dp);
    if (outerr == 0){
      std::cout << "end file output complete" << std::endl;
    }
    Kokkos::Profiling::popRegion();
  }
  Kokkos::finalize();
  return(0);
}

void periodic_cont(const int m, const int n, const int tl,
                   AllLevelType uAll, AllLevelType vAll, AllLevelType pAll){

  // North-South periodic continuation
  // Labeling conventions:
  // h = halo; i = interior;
  // north = m; south = 1;

  // Get desired time level of u, v, and p
  auto u = Kokkos::subview(uAll, Kokkos::ALL, tl);
  auto v = Kokkos::subview(vAll, Kokkos::ALL, tl);
  auto p = Kokkos::subview(pAll, Kokkos::ALL, tl);

  Kokkos::parallel_for("north-south_periodic_cont",
    Kokkos::RangePolicy<>(1,n+1),
    KOKKOS_LAMBDA(const int j){
      int hnorth = (m+1)*(n+2) + j;
      int hsouth = j;
      int inorth = m*(n+2) + j;
      int isouth = (n+2) + j;
      
      u(hsouth) = u(inorth);
      u(hnorth) = u(isouth);
      
      v(hsouth) = v(inorth);
      v(hnorth) = v(isouth);
      
      p(hsouth) = p(inorth);
      p(hnorth) = p(isouth);
    }
  );
  
  // East-West periodic continuation
  // h = halo; i = interior;
  // east = n; west = 1;
  
  // Corner point exchange indices
  
  int hsw = 0;
  int hne = (m+1)*(n+2)+n+1;
  int hse = n+1;
  int hnw = (m+1)*(n+2);
  
  int ine = m*(n+2)+n;
  int ise = 1*(n+2)+n;
  int isw = 1*(n+2)+1;
  int inw = m*(n+2)+1;
  
  Kokkos::parallel_for("east-west_periodic_cont",
    Kokkos::RangePolicy<>(1,m+1),
    KOKKOS_LAMBDA(const int i){
      int hwest = i*(n+2);
      int heast = i*(n+2) + n+1;
      int iwest = i*(n+2) + 1;
      int ieast = i*(n+2) + n;
      
      u(hwest) = u(ieast);
      u(heast) = u(iwest);
      
      v(hwest) = v(ieast);
      v(heast) = v(iwest);
      
      p(hwest) = p(ieast);
      p(heast) = p(iwest);
      
      if (i==m){
        u(hsw)   = u(ine);
        u(hnw)   = u(ise);
        u(hne)   = u(isw);
        u(hse)   = u(inw);

        v(hsw)   = v(ine);
        v(hnw)   = v(ise);
        v(hne)   = v(isw);
        v(hse)   = v(inw);
          
        p(hsw)   = p(ine);
        p(hnw)   = p(ise);
        p(hne)   = p(isw);
        p(hse)   = p(inw);
      }
    }
  );
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

// numerically the first step is special
struct first_step_functor{
  int m, n, tlmid, tlold, tlnew;
  AllLevelType uAll, vAll, pAll;
  double fsdx, fsdy, tdts8, tdtsdx, tdtsdy;

  first_step_functor(int _m, int _n, int _tlmid, int _tlold, int _tlnew,
                    AllLevelType _uAll, AllLevelType _vAll, AllLevelType _pAll,
                    double _fsdx, double _fsdy, double _tdts8, double _tdtsdx, double _tdtsdy){
    m      = _m;
    n      = _n;
    tlmid  = _tlmid;
    tlold  = _tlold;
    tlnew  = _tlnew;
    uAll   = _uAll;
    vAll   = _vAll;
    pAll   = _pAll;
    fsdx   = _fsdx;
    fsdy   = _fsdy;
    tdts8  = _tdts8;
    tdtsdx = _tdtsdx;
    tdtsdy = _tdtsdy;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const{
    auto u    = Kokkos::subview(uAll, Kokkos::ALL, tlmid);
    auto v    = Kokkos::subview(vAll, Kokkos::ALL, tlmid);
    auto p    = Kokkos::subview(pAll, Kokkos::ALL, tlmid);
    auto uold = Kokkos::subview(uAll, Kokkos::ALL, tlold);
    auto vold = Kokkos::subview(vAll, Kokkos::ALL, tlold);
    auto pold = Kokkos::subview(pAll, Kokkos::ALL, tlold);
    auto unew = Kokkos::subview(uAll, Kokkos::ALL, tlnew);
    auto vnew = Kokkos::subview(vAll, Kokkos::ALL, tlnew);
    auto pnew = Kokkos::subview(pAll, Kokkos::ALL, tlnew);

    int ij = i*(n+2)+j;
    int ijp1 = ij+1;
    int ijm1 = ij-1;
    
    int ip1j = ij+(n+2);
    int ip1jp1 = ip1j+1;
    int ip1jm1 = ip1j-1;
    
    int im1j= ij-(n+2);
    int im1jp1= im1j+1;
    
    double p00   = p(ij);      //10 mem access
    double pp0p1 = p(ijp1);    //6
    double pp1p0 = p(ip1j);    //6
    double pp1p1 = p(ip1jp1);  //3
    double u00   = u(ij);      //9
    double up0p1  = u(ijp1);   //4
    double um1p1 = u(im1jp1);  //4
    double um1p0 = u(im1j);    //4
    double v00   = v(ij);      //9
    double vp1p0 = v(ip1j);    //4
    double vp1m1 = v(ip1jm1);  //4
    double vp0m1 = v(ijm1);    //4
    
    double cut00 = .5 * (pp1p1 + pp0p1)       * up0p1;
    double cut10 = .5 * (pp0p1 + p(im1jp1)) * um1p1;
    double cut01 = .5 * (pp1p0 + p00)         * u00;
    double cut11 = .5 * (p00   + p(im1j))   * um1p0;
    
    double cvt00 = .5 * (pp1p1 + pp1p0)       * vp1p0;
    double cvt10 = .5 * (pp0p1 + p00)         * v00;
    double cvt01 = .5 * (pp1p0 + p(ip1jm1)) * vp1m1;
    double cvt11 = .5 * (p00   + p(ijm1))   * vp0m1;
    
    double zt00 = (fsdx * (vp1p0 - v00)     - fsdy * (up0p1 - u00))    /(p00     + pp1p0       + pp1p1 + pp0p1);
    double zt10 = (fsdx * (v00   - v(im1j)) - fsdy * (um1p1 - um1p0))  /(p(im1j) + p00         + pp0p1 + p(im1jp1));
    double zt01 = (fsdx * (vp1m1 - vp0m1)   - fsdy * (u00   - u(ijm1)))/(p(ijm1) + p(ip1jm1)   + pp1p0 + p00);
    
    double ht10 = pp0p1 + .25 * (up0p1     * up0p1     + um1p1 * um1p1 + v(ijp1) * v(ijp1)  + v00   * v00);
    double ht01 = pp1p0 + .25 * (u(ip1j) * u(ip1j) + u00   * u00   + vp1p0     * vp1p0      + vp1m1 * vp1m1);
    double ht11 = p00   + .25 * (u00       * u00       + um1p0 * um1p0 + v00       * v00        + vp0m1 * vp0m1);
    
    unew(ij) = uold(ij) + tdts8  * (zt00  + zt01) * (cvt00 + cvt10 + cvt11 + cvt01) - tdtsdx * (ht01 - ht11);
    vnew(ij) = vold(ij) - tdts8  * (zt00  + zt10) * (cut00 + cut10 + cut11 + cut01) - tdtsdy * (ht10 - ht11);
    pnew(ij) = pold(ij) - tdtsdx * (cut01 - cut11) - tdtsdy * (cvt10 - cvt11);
  }
};
void first_step(const int m, const int n, int tlmid, int tlold, int tlnew,
                AllLevelType uAll, AllLevelType vAll, AllLevelType pAll,
                double fsdx, double fsdy, double tdts8, double tdtsdx, double tdtsdy){

  first_step_functor functor(m, n, tlmid, tlold, tlnew,
                             uAll, vAll, pAll,
                             fsdx, fsdy, tdts8, tdtsdx, tdtsdy);

  Kokkos::parallel_for("first_step",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1,1},{m+1,n+1}),
    functor
  );
}

// update (assumes first step has been called).
struct update_functor{
  int m, n, tlmid, tlold, tlnew;
  AllLevelType uAll, vAll, pAll;
  double fsdx, fsdy, tdts8, tdtsdx, tdtsdy, alpha;

  update_functor(int _m, int _n, int _tlmid, int _tlold, int _tlnew,
                 AllLevelType _uAll, AllLevelType _vAll, AllLevelType _pAll,
                 double _fsdx, double _fsdy, double _tdts8, double _tdtsdx, double _tdtsdy, double _alpha){
    m      = _m;
    n      = _n;
    tlmid  = _tlmid;
    tlold  = _tlold;
    tlnew  = _tlnew;
    uAll   = _uAll;
    vAll   = _vAll;
    pAll   = _pAll;
    fsdx   = _fsdx;
    fsdy   = _fsdy;
    tdts8  = _tdts8;
    tdtsdx = _tdtsdx;
    tdtsdy = _tdtsdy;
    alpha  = _alpha;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const{
    auto u    = Kokkos::subview(uAll, Kokkos::ALL, tlmid);
    auto v    = Kokkos::subview(vAll, Kokkos::ALL, tlmid);
    auto p    = Kokkos::subview(pAll, Kokkos::ALL, tlmid);
    auto uold = Kokkos::subview(uAll, Kokkos::ALL, tlold);
    auto vold = Kokkos::subview(vAll, Kokkos::ALL, tlold);
    auto pold = Kokkos::subview(pAll, Kokkos::ALL, tlold);
    auto unew = Kokkos::subview(uAll, Kokkos::ALL, tlnew);
    auto vnew = Kokkos::subview(vAll, Kokkos::ALL, tlnew);
    auto pnew = Kokkos::subview(pAll, Kokkos::ALL, tlnew);

    int ij   = i*(n+2)+j;
    int ijp1 = ij+1;
    int ijm1 = ij-1;
    
    int ip1j   = ij+(n+2);
    int ip1jp1 = ip1j+1;
    int ip1jm1 = ip1j-1;
    
    int im1j   = ij-(n+2);
    int im1jp1 = im1j+1;
    
    double p00   = p(ij);      //10 mem access
    double pp0p1 = p(ijp1);    //6
    double pp1p0 = p(ip1j);    //6
    double pp1p1 = p(ip1jp1);  //3
    double u00   = u(ij);      //9
    double up0p1 = u(ijp1);    //4
    double um1p1 = u(im1jp1);  //4
    double um1p0 = u(im1j);    //4
    double v00   = v(ij);      //9
    double vp1p0 = v(ip1j);    //4
    double vp1m1 = v(ip1jm1);  //4
    double vp0m1 = v(ijm1);    //4
    
    double cut00 = .5 * (pp1p1 + pp0p1)     * up0p1;
    double cut10 = .5 * (pp0p1 + p(im1jp1)) * um1p1;
    double cut01 = .5 * (pp1p0 + p00)       * u00;
    double cut11 = .5 * (p00   + p(im1j))   * um1p0;
    
    double cvt00 = .5 * (pp1p1 + pp1p0)     * vp1p0;
    double cvt10 = .5 * (pp0p1 + p00)       * v00;
    double cvt01 = .5 * (pp1p0 + p(ip1jm1)) * vp1m1;
    double cvt11 = .5 * (p00   + p(ijm1))   * vp0m1;
    
    double zt00 = (fsdx * (vp1p0 - v00)     - fsdy * (up0p1 - u00))    /(p00     + pp1p0       + pp1p1 + pp0p1);
    double zt10 = (fsdx * (v00   - v(im1j)) - fsdy * (um1p1 - um1p0))  /(p(im1j) + p00         + pp0p1 + p(im1jp1));
    double zt01 = (fsdx * (vp1m1 - vp0m1)   - fsdy * (u00   - u(ijm1)))/(p(ijm1) + p(ip1jm1)   + pp1p0 + p00);
    
    double ht10 = pp0p1 + .25 * (up0p1     * up0p1     + um1p1 * um1p1 + v(ijp1) * v(ijp1)  + v00   * v00);
    double ht01 = pp1p0 + .25 * (u(ip1j) * u(ip1j) + u00   * u00   + vp1p0     * vp1p0      + vp1m1 * vp1m1);
    double ht11 = p00   + .25 * (u00       * u00       + um1p0 * um1p0 + v00       * v00        + vp0m1 * vp0m1);
    
    unew(ij) = uold(ij) + tdts8  * (zt00  + zt01) * (cvt00 + cvt10 + cvt11 + cvt01) - tdtsdx * (ht01 - ht11);
    vnew(ij) = vold(ij) - tdts8  * (zt00  + zt10) * (cut00 + cut10 + cut11 + cut01) - tdtsdy * (ht10 - ht11);
    pnew(ij) = pold(ij) - tdtsdx * (cut01 - cut11) - tdtsdy * (cvt10 - cvt11);
    
    uold(ij) = u00 + alpha * (unew(ij)  - 2. * u00 + uold(ij) );
    vold(ij) = v00 + alpha * (vnew(ij)  - 2. * v00 + vold(ij) );
    pold(ij) = p00 + alpha * (pnew(ij)  - 2. * p00 + pold(ij) );
  }
};
void update(const int m, const int n, int tlmid, int tlold, int tlnew,
            AllLevelType uAll, AllLevelType vAll, AllLevelType pAll,
            double fsdx, double fsdy, double tdts8, double tdtsdx, double tdtsdy, double alpha){

  update_functor functor(m, n, tlmid, tlold, tlnew,
                         uAll, vAll, pAll,
                         fsdx, fsdy, tdts8, tdtsdx, tdtsdy, alpha);
  Kokkos::parallel_for("update",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1,1},{m+1,n+1}),
    functor
  );
}
