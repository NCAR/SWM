// 1D range - Create i and j inside the flattened array
//   This updates all nodes except for the ghost ones
//   parallelization inside update()

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

void update(double u[DOMAIN_SIZE]) {
    
    auto R = range<1>{DOMAIN_SIZE};
    buffer<double, 1> u_buf(u, R);
    
    default_selector d_selector;
    queue q(d_selector);
    
    q.submit([&](handler &h) {
        
        auto u = u_buf.get_access(h, write_only);
        
        h.parallel_for(R, [=](auto index) { 
            int j = index%(N+2);
            int i = (int) (index - j)/(N+2);
            if (i==0 || j==0 || i == M+1 || j== N+1) {} 
            else {
                u[index] =  i*(N+2) + j;
            }
        }); 
    });
}

int main () {
    
    double u[DOMAIN_SIZE];
    for (int i=0;i<DOMAIN_SIZE; i++) u[i] = 0.;
    
    update(u);
    
    for (int i=0; i<DOMAIN_SIZE; i++)
        std::cout << "u[" << i << "] = " << u[i] << std::endl;
    
    return(0);
}

// Compile
// dpcpp -O2 -g -std=c++17 t005.cpp -o test005
