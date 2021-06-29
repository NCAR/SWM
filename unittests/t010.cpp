// Define buffers inside the main function and pass range and buffers to other functions

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
#define DIM 1

void update(queue q, range<DIM> R, buffer<double, DIM> u_buf) {
       
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
    
    default_selector d_selector;
    queue q(d_selector);
    std::cout << "Device: " << q.get_device().get_info<info::device::name>() << std::endl;
    
    double u[DOMAIN_SIZE];
    auto R = range<DIM>{DOMAIN_SIZE};
    buffer<double, DIM> u_buf(u, R);
    
    // Initialization
    for (int i=0;i<DOMAIN_SIZE; i++) u[i] = 0.;
    
    update(q, R, u_buf);
    
    // Synch
    host_accessor u_read(u_buf, read_only);
    
    for (int i=0; i<DOMAIN_SIZE; i++)
        std::cout << "u[" << i << "] = " << u[i] << std::endl;
    
    return(0);
}

// Compile
// dpcpp -O2 -g -std=c++17 t010.cpp -o test010
