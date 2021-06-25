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
            
            int j = ij%(N+2) + 1;
            int i = (int) (ij - j)/(N+2) + 1;
            
            int ipjp = (i+1)*(N+2) + j + 1;
            int ipj  = (i+1)*(N+2) + j;
            int ijp  = i*(N+2)     + j + 1;
            
            if (i >= M || j>= N) {} 
            else {
                u[ipjp] = -(psi[ipjp] - psi[ipj] - psi[ijp]) / (i*j);
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
    
    
    for (int i=0; i<DOMAIN_SIZE; i++)
        std::cout << "u[" << i << "] = " << u[i] << std::endl;
    
    return(0);
}

/* 
Compile
    dpcpp -O2 -g -std=c++17 t008.cpp -o test008
Run
    ./test008 > t008-check.sh
Check
    diff t008-check.sh t008.sh
*/