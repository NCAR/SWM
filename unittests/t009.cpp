// Test the indices for initial velocities

#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>
#include <cmath>
#include "dpc_common.hpp"
#include <limits>
#include <CL/sycl.hpp>
#if FPGA || FPGA_EMULATOR
#include <CL/sycl/INTEL/fpga_extensions.hpp>
#endif

using namespace std;
using namespace sycl;

#define M 4
#define N 5
#define DOMAIN_SIZE (M+2)*(N+2)

void update(queue q, double u[DOMAIN_SIZE], double psi[DOMAIN_SIZE]) {
    
    auto R = range<1>{DOMAIN_SIZE};
    buffer<double, 1> u_buf(u, R), psi_buf(psi, R);
    
    q.submit([&](handler &h) {
        
        auto u = u_buf.get_access(h, read_write);
        auto psi = psi_buf.get_access(h, write_only);
        
        h.parallel_for(R, [=](auto ij) {
            
            int j = ij%(N+2);
            int i = (int) (ij - j)/(N+2);
            
            int ijm1 = ij-1;
            int im1j= ij-(N+2);
            
            if (i==0 || j==0 || i == M+1 || j== N+1) {} 
            else {
                u[ij] = - .5 * (psi[ij] + 3.5 * psi[ijm1] - 5 * psi[im1j]) / (i*j);
            }
        }); 
    });
}

int main () {
    
    double u[DOMAIN_SIZE];
    for (int i=0;i<DOMAIN_SIZE; i++) u[i] = i;
    
    double psi[DOMAIN_SIZE];
    for (int i=0;i<DOMAIN_SIZE; i++) psi[i] = i*i;
    
    default_selector d_selector;
    queue q(d_selector);
    
    
    int n = N, m = M;
    update(q, u, psi);
    
    /* Test Created
    for (int ij=0; ij<DOMAIN_SIZE; ij++) {
        
        int j = ij%(n+2);
        int i = (int) (ij - j)/(n+2);
        
        int ijm1 = ij-1;
        int im1j= ij-(n+2);
        
        if (i==0 || j==0 || i == M+1 || j== N+1) {} 
        else {
            u[ij] = - .5 * (psi[ij] + 3.5 * psi[ijm1] - 5 * psi[im1j]) / (i*j);
        } 
    }
    
    // OR
    
    for (int i=1; i<m+1; i++) {
        for (int j=1;j<n+1;j++) {
          int ij = i*(n+2)+j;
          int ijm1 = ij-1;
          int im1j= ij-(n+2);
          
          u[ij] = - .5 * (psi[ij] + 3.5 * psi[ijm1] - 5 * psi[im1j]) / (i*j);
        }
    }
    */
    
    
    
    for (int i=0; i<DOMAIN_SIZE; i++)
        std::cout << "u[" << i << "] = " << u[i] << std::endl;
    
    return(0);
}

/* 
Compile
    dpcpp -O2 -g -std=c++17 t009.cpp -o test009
Run
    ./test009 > t009-check.sh
Check
    diff t009-check.sh t009.sh
*/